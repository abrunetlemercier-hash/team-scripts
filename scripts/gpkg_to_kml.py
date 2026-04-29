#!/usr/bin/env python3
"""
Convert GeoPackage feature tables to KML.

The converter is intentionally dependency-free so it can run in the same
lightweight app environment as the KML to GeoPackage tool.
"""

import io
import os
import shutil
import sqlite3
import struct
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from scripts._base import ScriptResult

KML_NS = "http://www.opengis.net/kml/2.2"
ET.register_namespace("", KML_NS)


def _read_uint32(blob: bytes, offset: int, endian: str) -> tuple[int, int]:
    return struct.unpack_from(f"{endian}I", blob, offset)[0], offset + 4


def _read_double(blob: bytes, offset: int, endian: str) -> tuple[float, int]:
    return struct.unpack_from(f"{endian}d", blob, offset)[0], offset + 8


def _gpkg_wkb_offset(blob: bytes) -> int:
    """Return the WKB start offset inside a GeoPackage geometry blob."""
    if len(blob) < 8 or blob[:2] != b"GP":
        return 0

    flags = blob[3]
    envelope_type = (flags >> 1) & 0b111
    envelope_sizes = {
        0: 0,
        1: 32,
        2: 48,
        3: 48,
        4: 64,
    }
    return 8 + envelope_sizes.get(envelope_type, 0)


def _parse_wkb_geometry(blob: bytes):
    """Parse WKB Polygon/MultiPolygon geometry into nested coordinate rings."""
    offset = _gpkg_wkb_offset(blob)
    if offset >= len(blob):
        return None

    byte_order = blob[offset]
    endian = "<" if byte_order == 1 else ">"
    offset += 1
    geom_type, offset = _read_uint32(blob, offset, endian)
    geom_type = geom_type % 1000

    def read_point(pos: int) -> tuple[tuple[float, float], int]:
        lon, pos = _read_double(blob, pos, endian)
        lat, pos = _read_double(blob, pos, endian)
        return (lon, lat), pos

    def read_polygon(pos: int) -> tuple[list[list[tuple[float, float]]], int]:
        ring_count, pos = _read_uint32(blob, pos, endian)
        rings = []
        for _ in range(ring_count):
            point_count, pos = _read_uint32(blob, pos, endian)
            ring = []
            for _ in range(point_count):
                point, pos = read_point(pos)
                ring.append(point)
            rings.append(ring)
        return rings, pos

    if geom_type == 3:
        rings, _ = read_polygon(offset)
        return [rings]

    if geom_type == 6:
        polygon_count, offset = _read_uint32(blob, offset, endian)
        polygons = []
        for _ in range(polygon_count):
            nested_order = blob[offset]
            nested_endian = "<" if nested_order == 1 else ">"
            if nested_endian != endian:
                endian = nested_endian
            offset += 1
            nested_type, offset = _read_uint32(blob, offset, endian)
            if nested_type % 1000 != 3:
                return None
            rings, offset = read_polygon(offset)
            polygons.append(rings)
        return polygons

    return None


def _coordinates_text(ring: list[tuple[float, float]]) -> str:
    return " ".join(f"{lon:.8f},{lat:.8f},0" for lon, lat in ring)


def _add_polygon(parent, polygons):
    if len(polygons) > 1:
        multi = ET.SubElement(parent, f"{{{KML_NS}}}MultiGeometry")
        for rings in polygons:
            polygon = ET.SubElement(multi, f"{{{KML_NS}}}Polygon")
            _add_polygon_rings(polygon, rings)
        return

    polygon = ET.SubElement(parent, f"{{{KML_NS}}}Polygon")
    _add_polygon_rings(polygon, polygons[0])


def _add_polygon_rings(polygon, rings):
    if not rings:
        return

    outer = ET.SubElement(polygon, f"{{{KML_NS}}}outerBoundaryIs")
    outer_ring = ET.SubElement(outer, f"{{{KML_NS}}}LinearRing")
    outer_coords = ET.SubElement(outer_ring, f"{{{KML_NS}}}coordinates")
    outer_coords.text = _coordinates_text(rings[0])

    for ring in rings[1:]:
        inner = ET.SubElement(polygon, f"{{{KML_NS}}}innerBoundaryIs")
        inner_ring = ET.SubElement(inner, f"{{{KML_NS}}}LinearRing")
        inner_coords = ET.SubElement(inner_ring, f"{{{KML_NS}}}coordinates")
        inner_coords.text = _coordinates_text(ring)


def _safe_layer_name(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in (" ", "-", "_") else "_" for ch in name).strip()


def _feature_tables(conn: sqlite3.Connection) -> list[str]:
    rows = conn.execute(
        "SELECT table_name FROM gpkg_contents WHERE data_type = 'features' ORDER BY table_name"
    ).fetchall()
    if rows:
        return [r[0] for r in rows]

    return [
        r[0]
        for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'gpkg_%' "
            "AND name NOT LIKE 'rtree_%' ORDER BY name"
        ).fetchall()
    ]


def convert_gpkg_to_kml(gpkg_path: Path, output_dir: Path) -> tuple[int, list[Path], str]:
    """Convert one GPKG file. Returns (feature_count, output_paths, log)."""
    conn = sqlite3.connect(str(gpkg_path))
    conn.row_factory = sqlite3.Row
    output_paths = []
    logs = []
    total_features = 0

    for table in _feature_tables(conn):
        columns = conn.execute(f'PRAGMA table_info("{table}")').fetchall()
        column_names = [c["name"] for c in columns]
        geom_col = "geom" if "geom" in column_names else None
        if geom_col is None:
            for col in column_names:
                if col.lower() in {"geometry", "wkb_geometry"}:
                    geom_col = col
                    break
        if geom_col is None:
            logs.append(f"  {table}: skipped, no geometry column")
            continue

        rows = conn.execute(f'SELECT * FROM "{table}"').fetchall()
        if not rows:
            logs.append(f"  {table}: skipped, no features")
            continue

        kml = ET.Element(f"{{{KML_NS}}}kml")
        document = ET.SubElement(kml, f"{{{KML_NS}}}Document")
        doc_name = ET.SubElement(document, f"{{{KML_NS}}}name")
        layer_name = _safe_layer_name(table) or "features"
        doc_name.text = f"{gpkg_path.stem} - {layer_name}"

        written = 0
        for idx, row in enumerate(rows, start=1):
            geom_blob = row[geom_col]
            if geom_blob is None:
                continue
            polygons = _parse_wkb_geometry(bytes(geom_blob))
            if not polygons:
                continue

            placemark = ET.SubElement(document, f"{{{KML_NS}}}Placemark")
            row_name = row["Name"] if "Name" in row.keys() and row["Name"] else f"{table} {idx}"
            name_el = ET.SubElement(placemark, f"{{{KML_NS}}}name")
            name_el.text = str(row_name)

            extended = ET.SubElement(placemark, f"{{{KML_NS}}}ExtendedData")
            for col in column_names:
                if col in {geom_col, "fid"}:
                    continue
                value = row[col]
                if value is None:
                    continue
                data = ET.SubElement(extended, f"{{{KML_NS}}}Data", name=col)
                value_el = ET.SubElement(data, f"{{{KML_NS}}}value")
                value_el.text = str(value)

            _add_polygon(placemark, polygons)
            written += 1

        if written == 0:
            logs.append(f"  {table}: skipped, no supported Polygon/MultiPolygon features")
            continue

        suffix = "" if len(_feature_tables(conn)) == 1 else f"_{layer_name}"
        kml_path = output_dir / f"{gpkg_path.stem}{suffix}.kml"
        ET.indent(kml, space="  ")
        ET.ElementTree(kml).write(kml_path, encoding="utf-8", xml_declaration=True)
        output_paths.append(kml_path)
        total_features += written
        logs.append(f"  {table}: {written} features -> {kml_path.name}")

    conn.close()
    return total_features, output_paths, "\n".join(logs)


def run_conversion(base_dir: Path) -> tuple[str, str | None]:
    """Run GPKG to KML conversion. Returns (log output, zip path)."""
    buf = io.StringIO()
    gpkg_files = sorted(base_dir.glob("*.gpkg"))

    if not gpkg_files:
        return "No GPKG files found. Upload .gpkg files first.", None

    output_dir = base_dir / "kml_output"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    old_zip = base_dir / "kml_files.zip"
    if old_zip.exists():
        old_zip.unlink()

    buf.write(f"Found {len(gpkg_files)} GPKG files\n\n")
    all_outputs = []
    total = 0
    for gpkg_path in gpkg_files:
        buf.write(f"Processing: {gpkg_path.name}\n")
        try:
            count, paths, log = convert_gpkg_to_kml(gpkg_path, output_dir)
            total += count
            all_outputs.extend(paths)
            if log:
                buf.write(log + "\n")
            buf.write("\n")
        except Exception as exc:
            buf.write(f"  Error: {exc}\n\n")

    if not all_outputs:
        return buf.getvalue() + "No KML files were created.", None

    zip_path = base_dir / "kml_files.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in all_outputs:
            zf.write(path, path.name)

    buf.write("=" * 60 + "\n")
    buf.write(f"Done! {total} features converted.\n")
    buf.write(f"{len(all_outputs)} KML file(s) ready for download.")
    return buf.getvalue(), str(zip_path)


class GpkgToKmlScript:
    name = "GPKG to KML"
    description = "Convert uploaded GeoPackage files back to KML while preserving feature attributes."
    input_type = "gpkg"

    async def run(self) -> ScriptResult:
        try:
            from app.config import GPKG_DIR

            GPKG_DIR.mkdir(parents=True, exist_ok=True)
            output, zip_path = run_conversion(GPKG_DIR)

            for gpkg_file in GPKG_DIR.glob("*.gpkg"):
                gpkg_file.unlink()
            output_dir = GPKG_DIR / "kml_output"
            if output_dir.exists():
                shutil.rmtree(output_dir)

            return ScriptResult(
                success=zip_path is not None,
                output=output,
                error=None if zip_path else "No KML files were created.",
                download_file=zip_path,
            )
        except Exception as exc:
            return ScriptResult(success=False, output="", error=str(exc))


script = GpkgToKmlScript()


if __name__ == "__main__":
    from app.config import GPKG_DIR

    GPKG_DIR.mkdir(parents=True, exist_ok=True)
    output, _ = run_conversion(GPKG_DIR)
    print(output)
