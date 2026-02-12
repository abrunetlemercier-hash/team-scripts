#!/usr/bin/env python3
"""
Convert KML files to GeoPackage with clean attribute extraction.

Handles 3 KML variants:
  1. ExtendedData/SimpleData (structured attributes)
  2. HTML description tables (attributes in <td> pairs)
  3. Geometry-only (no attributes, just name)

Post-processing for the 3 main Java files (East/Central/West):
  - Rename polygons (Name) → Province name
  - Status: "Private" if original name contained "Private", else "APL"
  - Land_Type: extract from original name if hint present (Pond/Coastal),
    otherwise keep source value or NULL (no extrapolation)
  - Province: fill from file context if NULL
  - District, Sub_Distri, Village: fill NULLs by nearest neighbor (haversine)
  - Fix typos: "Ponds" → "Pond"

Geometry fields (always computed from geometry, never from description):
  - Longitude / Latitude: centroid of polygon in EPSG:4326
  - Area_Ha: polygon area reprojected to UTM 48S (EPSG:32748), in hectares

Output:
  - One .gpkg per KML (Java files renamed by province)
  - ALL_merged.gpkg with only the 3 Java files + Source_File column
"""

import io
import xml.etree.ElementTree as ET
import re
import sqlite3
import struct
import math
import os
import sys
from pathlib import Path

from scripts._base import ScriptResult

# ─── Configuration ───────────────────────────────────────────────────────────

# Target attribute fields (text)
TEXT_FIELDS = [
    "Province", "District", "Sub_Distri", "Village",
    "Land_Type", "Status", "Code", "NRegen_Typ", "ARR_type", "Accretion",
]

# Extra text fields found in some files
EXTRA_TEXT_FIELDS = [
    "M_Species", "P_Method", "N_Seeds", "Group_Name", "Group_Lead", "PLOT_ID",
    "Rev_1212", "Source_File",
]

# All text fields we look for
ALL_TEXT_FIELDS = TEXT_FIELDS + EXTRA_TEXT_FIELDS

# Numeric fields computed from geometry
GEOM_FIELDS = ["Longitude", "Latitude", "Area_Ha"]

KML_NS = "{http://www.opengis.net/kml/2.2}"

# Map from KML stem to clean province-based GPKG name
RENAME_MAP = {
    "20260126_East_Java": "East Java",
    "20260126_Central_Java": "Central Java",
    "20260126_West_Java": "West Java",
}

# KML stems that get the full Java post-processing
JAVA_STEMS = set(RENAME_MAP.keys())

# ─── Geometry helpers ────────────────────────────────────────────────────────

def parse_coordinates(coord_text):
    """Parse KML coordinate string into list of (lon, lat) tuples."""
    coords = []
    for token in coord_text.strip().split():
        parts = token.split(",")
        if len(parts) >= 2:
            lon, lat = float(parts[0]), float(parts[1])
            coords.append((lon, lat))
    return coords


def polygon_centroid(rings):
    """Compute centroid of polygon (first ring = exterior).
    Uses the signed-area formula for 2D polygons."""
    if not rings or not rings[0]:
        return (0, 0)
    ring = rings[0]
    n = len(ring)
    if n < 3:
        lons = [p[0] for p in ring]
        lats = [p[1] for p in ring]
        return (sum(lons) / len(lons), sum(lats) / len(lats))

    area = 0.0
    cx = 0.0
    cy = 0.0
    for i in range(n - 1):
        x0, y0 = ring[i]
        x1, y1 = ring[i + 1]
        cross = x0 * y1 - x1 * y0
        area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross

    area *= 0.5
    if abs(area) < 1e-15:
        lons = [p[0] for p in ring]
        lats = [p[1] for p in ring]
        return (sum(lons) / len(lons), sum(lats) / len(lats))

    cx /= (6.0 * area)
    cy /= (6.0 * area)
    return (cx, cy)


def latlon_to_utm48s(lon, lat):
    """Convert WGS84 (lon, lat) to UTM zone 48S (EPSG:32748).
    Returns (easting, northing)."""
    lon0 = 105.0  # central meridian for UTM zone 48
    k0 = 0.9996
    a = 6378137.0
    f = 1 / 298.257223563
    e2 = 2 * f - f * f
    e_prime2 = e2 / (1 - e2)

    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    lon0_rad = math.radians(lon0)

    N = a / math.sqrt(1 - e2 * math.sin(lat_rad) ** 2)
    T = math.tan(lat_rad) ** 2
    C = e_prime2 * math.cos(lat_rad) ** 2
    A = (lon_rad - lon0_rad) * math.cos(lat_rad)

    e4 = e2 * e2
    e6 = e4 * e2
    M = a * (
        (1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * lat_rad
        - (3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * math.sin(2 * lat_rad)
        + (15 * e4 / 256 + 45 * e6 / 1024) * math.sin(4 * lat_rad)
        - (35 * e6 / 3072) * math.sin(6 * lat_rad)
    )

    A2 = A * A
    A4 = A2 * A2
    A6 = A4 * A2

    easting = k0 * N * (
        A + (1 - T + C) * A2 * A / 6
        + (5 - 18 * T + T * T + 72 * C - 58 * e_prime2) * A4 * A / 120
    ) + 500000.0

    northing = k0 * (
        M + N * math.tan(lat_rad) * (
            A2 / 2
            + (5 - T + 9 * C + 4 * C * C) * A4 / 24
            + (61 - 58 * T + T * T + 600 * C - 330 * e_prime2) * A6 / 720
        )
    )

    if lat < 0:
        northing += 10000000.0

    return (easting, northing)


def polygon_area_utm(rings):
    """Compute polygon area in m² by projecting all rings to UTM 48S."""
    total = 0.0
    for idx, ring in enumerate(rings):
        projected = [latlon_to_utm48s(lon, lat) for lon, lat in ring]
        n = len(projected)
        area = 0.0
        for i in range(n - 1):
            x0, y0 = projected[i]
            x1, y1 = projected[i + 1]
            area += x0 * y1 - x1 * y0
        area = abs(area) / 2.0
        if idx == 0:
            total += area  # exterior
        else:
            total -= area  # holes
    return abs(total)


def rings_to_gpkg_wkb(rings):
    """Encode polygon rings as GeoPackage WKB."""
    wkb = bytearray()
    wkb.append(0x01)  # little-endian
    wkb.extend(struct.pack("<I", 3))  # wkbPolygon

    wkb.extend(struct.pack("<I", len(rings)))
    for ring in rings:
        wkb.extend(struct.pack("<I", len(ring)))
        for lon, lat in ring:
            wkb.extend(struct.pack("<dd", lon, lat))

    all_coords = [c for ring in rings for c in ring]
    lons = [c[0] for c in all_coords]
    lats = [c[1] for c in all_coords]
    minx, maxx = min(lons), max(lons)
    miny, maxy = min(lats), max(lats)

    header = bytearray()
    header.extend(b"GP")
    header.append(0x00)
    header.append(0b00000011)  # LE + envelope type 1
    header.extend(struct.pack("<i", 4326))
    header.extend(struct.pack("<dddd", minx, maxx, miny, maxy))

    return bytes(header) + bytes(wkb)


def multipolygon_to_gpkg_wkb(polygons_rings):
    """Encode multipolygon as GeoPackage WKB."""
    wkb = bytearray()
    wkb.append(0x01)
    wkb.extend(struct.pack("<I", 6))  # wkbMultiPolygon
    wkb.extend(struct.pack("<I", len(polygons_rings)))

    for rings in polygons_rings:
        wkb.append(0x01)
        wkb.extend(struct.pack("<I", 3))
        wkb.extend(struct.pack("<I", len(rings)))
        for ring in rings:
            wkb.extend(struct.pack("<I", len(ring)))
            for lon, lat in ring:
                wkb.extend(struct.pack("<dd", lon, lat))

    all_coords = [c for pr in polygons_rings for ring in pr for c in ring]
    lons = [c[0] for c in all_coords]
    lats = [c[1] for c in all_coords]

    header = bytearray()
    header.extend(b"GP")
    header.append(0x00)
    header.append(0b00000011)
    header.extend(struct.pack("<i", 4326))
    header.extend(struct.pack("<dddd", min(lons), max(lons), min(lats), max(lats)))

    return bytes(header) + bytes(wkb)


# ─── KML Parsing ─────────────────────────────────────────────────────────────

def extract_simpledata(placemark):
    """Extract attributes from ExtendedData/SchemaData/SimpleData."""
    attrs = {}
    for sd in placemark.iter(f"{KML_NS}SimpleData"):
        name = sd.get("name")
        val = (sd.text or "").strip()
        if name:
            attrs[name] = val if val else None
    return attrs


def extract_html_description(placemark):
    """Extract attributes from HTML description table (<td>Key</td><td>Value</td>)."""
    attrs = {}
    desc_el = placemark.find(f"{KML_NS}description")
    if desc_el is None or not desc_el.text:
        return attrs

    desc = desc_el.text
    for field in ALL_TEXT_FIELDS + ["FID", "fid", "Area_Ha", "Longitude", "Latitude"]:
        pattern = (
            r'(?si)<td[^>]*>\s*(?:<b>)?\s*'
            + re.escape(field)
            + r'\s*(?:</b>)?\s*</td>\s*<td[^>]*>(.*?)</td>'
        )
        m = re.search(pattern, desc)
        if m:
            val = m.group(1).strip()
            attrs[field] = val if val else None

    return attrs


def extract_polygon_rings(polygon_el):
    """Extract rings from a KML Polygon element."""
    rings = []
    outer = polygon_el.find(f"{KML_NS}outerBoundaryIs")
    if outer is not None:
        lr = outer.find(f"{KML_NS}LinearRing")
        if lr is not None:
            coord_el = lr.find(f"{KML_NS}coordinates")
            if coord_el is not None and coord_el.text:
                rings.append(parse_coordinates(coord_el.text))

    for inner in polygon_el.findall(f"{KML_NS}innerBoundaryIs"):
        lr = inner.find(f"{KML_NS}LinearRing")
        if lr is not None:
            coord_el = lr.find(f"{KML_NS}coordinates")
            if coord_el is not None and coord_el.text:
                rings.append(parse_coordinates(coord_el.text))

    return rings


def parse_placemark(pm):
    """Parse a single Placemark. Returns (attrs_dict, polygon_rings_list)."""
    attrs = extract_simpledata(pm)
    if not attrs:
        attrs = extract_html_description(pm)

    name_el = pm.find(f"{KML_NS}name")
    if name_el is not None and name_el.text:
        attrs.setdefault("Name", name_el.text.strip())

    all_polygon_rings = []
    polygon_el = pm.find(f"{KML_NS}Polygon")
    if polygon_el is not None:
        rings = extract_polygon_rings(polygon_el)
        if rings:
            all_polygon_rings.append(rings)

    mg = pm.find(f"{KML_NS}MultiGeometry")
    if mg is not None:
        for poly in mg.findall(f"{KML_NS}Polygon"):
            rings = extract_polygon_rings(poly)
            if rings:
                all_polygon_rings.append(rings)

    return attrs, all_polygon_rings


def parse_kml(filepath):
    """Parse a KML file, return list of (attrs_dict, polygon_rings_list)."""
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
    # Fix missing xsi namespace declaration
    if "xsi:" in content and 'xmlns:xsi' not in content:
        content = content.replace(
            '<kml ',
            '<kml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
            1
        )
    root = ET.fromstring(content)

    features = []
    for pm in root.iter(f"{KML_NS}Placemark"):
        attrs, polygon_rings = parse_placemark(pm)
        if polygon_rings:
            features.append((attrs, polygon_rings))

    return features


# ─── GeoPackage Writer ───────────────────────────────────────────────────────

def create_gpkg(filepath, features, source_name):
    """Create a GeoPackage from parsed features."""
    if os.path.exists(filepath):
        os.remove(filepath)

    conn = sqlite3.connect(filepath)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA application_id=0x47504B47")
    conn.execute("PRAGMA user_version=10301")

    # Determine which text fields have data
    fields_with_data = set()
    for attrs, _ in features:
        for key in attrs:
            if key not in ("FID", "fid", "Area_Ha", "Longitude", "Latitude", "Name"):
                fields_with_data.add(key)

    text_fields = [f for f in ALL_TEXT_FIELDS if f in fields_with_data]
    has_name = any("Name" in attrs for attrs, _ in features)
    if has_name:
        text_fields.insert(0, "Name")

    # ── GeoPackage metadata tables ──
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS gpkg_spatial_ref_sys (
            srs_name TEXT NOT NULL, srs_id INTEGER NOT NULL PRIMARY KEY,
            organization TEXT NOT NULL, organization_coordsys_id INTEGER NOT NULL,
            definition TEXT NOT NULL, description TEXT);

        CREATE TABLE IF NOT EXISTS gpkg_contents (
            table_name TEXT NOT NULL PRIMARY KEY,
            data_type TEXT NOT NULL DEFAULT 'features',
            identifier TEXT UNIQUE, description TEXT DEFAULT '',
            last_change DATETIME NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')),
            min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE,
            srs_id INTEGER,
            CONSTRAINT fk_gc_r_srs_id FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys(srs_id));

        CREATE TABLE IF NOT EXISTS gpkg_geometry_columns (
            table_name TEXT NOT NULL, column_name TEXT NOT NULL,
            geometry_type_name TEXT NOT NULL, srs_id INTEGER NOT NULL,
            z TINYINT NOT NULL, m TINYINT NOT NULL,
            CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name),
            CONSTRAINT fk_gc_tn FOREIGN KEY (table_name) REFERENCES gpkg_contents(table_name),
            CONSTRAINT fk_gc_srs FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys(srs_id));
    """)

    conn.execute("""
        INSERT OR IGNORE INTO gpkg_spatial_ref_sys
        (srs_name, srs_id, organization, organization_coordsys_id, definition)
        VALUES (?, ?, ?, ?, ?)
    """, (
        "WGS 84", 4326, "EPSG", 4326,
        'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],'
        'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],'
        'AUTHORITY["EPSG","4326"]]'
    ))
    conn.execute("""INSERT OR IGNORE INTO gpkg_spatial_ref_sys
        (srs_name, srs_id, organization, organization_coordsys_id, definition)
        VALUES ('Undefined cartesian SRS', -1, 'NONE', -1, 'undefined')""")
    conn.execute("""INSERT OR IGNORE INTO gpkg_spatial_ref_sys
        (srs_name, srs_id, organization, organization_coordsys_id, definition)
        VALUES ('Undefined geographic SRS', 0, 'NONE', 0, 'undefined')""")

    # ── Feature table ──
    tbl = "features"
    col_defs = ["fid INTEGER PRIMARY KEY AUTOINCREMENT", "geom BLOB"]
    for f in text_fields:
        col_defs.append(f'"{f}" TEXT')
    col_defs.append("Longitude REAL")
    col_defs.append("Latitude REAL")
    col_defs.append("Area_Ha REAL")
    conn.execute(f'CREATE TABLE "{tbl}" ({", ".join(col_defs)})')

    # Envelope
    all_lons, all_lats = [], []
    for _, polygon_rings in features:
        for rings in polygon_rings:
            for ring in rings:
                for lon, lat in ring:
                    all_lons.append(lon)
                    all_lats.append(lat)
    minx = min(all_lons) if all_lons else 0
    maxx = max(all_lons) if all_lons else 0
    miny = min(all_lats) if all_lats else 0
    maxy = max(all_lats) if all_lats else 0

    has_multi = any(len(pr) > 1 for _, pr in features)
    geom_type = "MULTIPOLYGON" if has_multi else "POLYGON"

    conn.execute("""INSERT INTO gpkg_contents
        (table_name, data_type, identifier, description, min_x, min_y, max_x, max_y, srs_id)
        VALUES (?, 'features', ?, ?, ?, ?, ?, ?, 4326)
    """, (tbl, source_name, f"Converted from {source_name}", minx, miny, maxx, maxy))

    conn.execute("""INSERT INTO gpkg_geometry_columns
        (table_name, column_name, geometry_type_name, srs_id, z, m)
        VALUES (?, 'geom', ?, 4326, 0, 0)
    """, (tbl, geom_type))

    # ── Insert features ──
    col_names = ["geom"] + [f'"{f}"' for f in text_fields] + ["Longitude", "Latitude", "Area_Ha"]
    placeholders = ", ".join(["?"] * len(col_names))
    insert_sql = f'INSERT INTO "{tbl}" ({", ".join(col_names)}) VALUES ({placeholders})'

    rows = []
    for attrs, polygon_rings in features:
        if len(polygon_rings) == 1 and not has_multi:
            wkb = rings_to_gpkg_wkb(polygon_rings[0])
        else:
            wkb = multipolygon_to_gpkg_wkb(polygon_rings)

        all_exterior = [rings[0] for rings in polygon_rings if rings]
        if len(all_exterior) == 1:
            cx, cy = polygon_centroid([all_exterior[0]])
        else:
            total_a, cx, cy = 0, 0, 0
            for ext in all_exterior:
                c = polygon_centroid([ext])
                a = polygon_area_utm([ext])
                cx += c[0] * a
                cy += c[1] * a
                total_a += a
            if total_a > 0:
                cx /= total_a
                cy /= total_a

        longitude = round(cx, 6)
        latitude = round(cy, 6)
        area_ha = round(sum(polygon_area_utm(rings) for rings in polygon_rings) / 10000.0, 4)

        text_vals = []
        for f in text_fields:
            val = attrs.get(f)
            if val is not None:
                val = val.strip()
                if not val:
                    val = None
            text_vals.append(val)

        rows.append([wkb] + text_vals + [longitude, latitude, area_ha])

    conn.executemany(insert_sql, rows)
    conn.commit()

    # Spatial index
    conn.execute("""CREATE TABLE "rtree_features_geom" (
        id INTEGER PRIMARY KEY, minx REAL, maxx REAL, miny REAL, maxy REAL)""")
    for row in conn.execute(f'SELECT fid, geom FROM "{tbl}"'):
        fid, geom_blob = row
        if geom_blob and len(geom_blob) > 40:
            env = struct.unpack_from("<dddd", geom_blob, 8)
            conn.execute(
                'INSERT INTO "rtree_features_geom" VALUES (?,?,?,?,?)',
                (fid, env[0], env[1], env[2], env[3]))
    conn.commit()
    conn.close()
    return len(rows), text_fields


# ─── Post-processing: Java files ─────────────────────────────────────────────

def haversine(lon1, lat1, lon2, lat2):
    """Distance in km between two WGS84 points."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return R * 2 * math.asin(math.sqrt(a))


def fill_nulls_by_nearest(conn, field):
    """Fill NULL values of `field` by copying from nearest feature that has it."""
    ref_rows = conn.execute(
        f'SELECT fid, Longitude, Latitude, "{field}" '
        f'FROM features WHERE "{field}" IS NOT NULL'
    ).fetchall()
    null_rows = conn.execute(
        f'SELECT fid, Longitude, Latitude '
        f'FROM features WHERE "{field}" IS NULL'
    ).fetchall()
    count = 0
    for fid, lon, lat in null_rows:
        best_dist = float('inf')
        best_val = None
        for _, ref_lon, ref_lat, ref_val in ref_rows:
            d = haversine(lon, lat, ref_lon, ref_lat)
            if d < best_dist:
                best_dist = d
                best_val = ref_val
        if best_val is not None:
            conn.execute(f'UPDATE features SET "{field}" = ? WHERE fid = ?',
                         (best_val, fid))
            count += 1
    return count


def postprocess_java_file(gpkg_path, province, kml_path):
    """Apply all post-processing rules to a Java province GPKG.

    Rules:
      1. Extract Land_Type from original Name if hint present (Pond/Coastal),
         but do NOT extrapolate by proximity — leave NULL if unknown.
      2. Fix typo "Ponds" → "Pond".
      3. Status: "Private" if original name contains "private", else "APL".
      4. Rename all polygon Names → province name.
      5. Fill NULL Province with province name.
      6. Fill NULL District, Sub_Distri, Village by nearest neighbor.
    """
    # Re-parse KML to get original Name and Land_Type per feature (in order)
    features = parse_kml(str(kml_path))
    originals = []
    for attrs, _ in features:
        orig_name = attrs.get("Name", "") or ""
        orig_lt = (attrs.get("Land_Type") or "").strip() or None
        originals.append((orig_name, orig_lt))

    conn = sqlite3.connect(str(gpkg_path))

    # ── 1. Land_Type: extract from original name where source was NULL ──
    for i, (orig_name, orig_lt) in enumerate(originals):
        fid = i + 1
        name_lower = orig_name.lower()
        if orig_lt is None:
            # Check name for hints
            if "pond" in name_lower:
                conn.execute('UPDATE features SET Land_Type = ? WHERE fid = ?',
                             ("Pond", fid))
            elif "coastal" in name_lower:
                conn.execute('UPDATE features SET Land_Type = ? WHERE fid = ?',
                             ("Coastal", fid))
            # else: leave NULL — no extrapolation

    # ── 2. Fix typo ──
    conn.execute("UPDATE features SET Land_Type = 'Pond' WHERE Land_Type = 'Ponds'")

    # ── 3. Status from original name ──
    for i, (orig_name, _) in enumerate(originals):
        fid = i + 1
        if "private" in orig_name.lower():
            conn.execute("UPDATE features SET Status = 'Private' WHERE fid = ?", (fid,))
        else:
            conn.execute("UPDATE features SET Status = 'APL' WHERE fid = ?", (fid,))

    # ── 4. Rename all polygon Names → province ──
    conn.execute("UPDATE features SET Name = ?", (province,))

    # ── 5. Fill NULL Province ──
    conn.execute("UPDATE features SET Province = ? WHERE Province IS NULL", (province,))

    # ── 6. Fill NULL District, Sub_Distri, Village by nearest neighbor ──
    conn.commit()
    for field in ["District", "Sub_Distri", "Village"]:
        n = fill_nulls_by_nearest(conn, field)
        if n > 0:
            print(f"    {field}: {n} NULLs filled by nearest neighbor")

    conn.commit()

    # ── Report ──
    total = conn.execute("SELECT COUNT(*) FROM features").fetchone()[0]
    prv = conn.execute("SELECT COUNT(*) FROM features WHERE Status='Private'").fetchone()[0]
    apl = conn.execute("SELECT COUNT(*) FROM features WHERE Status='APL'").fetchone()[0]
    lt_null = conn.execute("SELECT COUNT(*) FROM features WHERE Land_Type IS NULL").fetchone()[0]
    nulls = {}
    for col in ["Province", "District", "Sub_Distri", "Village"]:
        n = conn.execute(f'SELECT COUNT(*) FROM features WHERE "{col}" IS NULL').fetchone()[0]
        if n > 0:
            nulls[col] = n
    print(f"    {total} features | Private={prv}, APL={apl} | Land_Type NULL={lt_null}")
    if nulls:
        print(f"    Remaining NULLs: {nulls}")

    conn.close()


# ─── Merged GPKG builder ─────────────────────────────────────────────────────

def create_merged_gpkg(output_path, source_gpkgs):
    """Create ALL_merged.gpkg from a list of (gpkg_path, source_label) tuples."""
    # Get union of all columns
    all_col_info = {}
    for gpkg_path, _ in source_gpkgs:
        conn = sqlite3.connect(str(gpkg_path))
        for c in conn.execute('PRAGMA table_info(features)').fetchall():
            _, name, ctype, _, _, _ = c
            if name != 'fid':
                all_col_info[name] = ctype
        conn.close()

    non_fid_cols = list(all_col_info.items())
    col_names = [n for n, _ in non_fid_cols]

    # Read all features with unified columns
    all_rows = []
    for gpkg_path, src_label in source_gpkgs:
        conn = sqlite3.connect(str(gpkg_path))
        src_cols = [c[1] for c in conn.execute('PRAGMA table_info(features)').fetchall()]
        select_parts = []
        for name in col_names:
            if name in src_cols:
                select_parts.append(f'"{name}"')
            else:
                select_parts.append(f'NULL as "{name}"')
        rows = conn.execute(f'SELECT {", ".join(select_parts)} FROM features').fetchall()
        for r in rows:
            all_rows.append(list(r) + [src_label])
        conn.close()

    col_names.append('Source_File')
    non_fid_cols.append(('Source_File', 'TEXT'))

    if os.path.exists(str(output_path)):
        os.remove(str(output_path))

    conn = sqlite3.connect(str(output_path))
    conn.execute('PRAGMA journal_mode=WAL')
    conn.execute('PRAGMA application_id=0x47504B47')
    conn.execute('PRAGMA user_version=10301')

    conn.execute("""CREATE TABLE gpkg_spatial_ref_sys (
        srs_name TEXT NOT NULL, srs_id INTEGER NOT NULL PRIMARY KEY,
        organization TEXT NOT NULL, organization_coordsys_id INTEGER NOT NULL,
        definition TEXT NOT NULL, description TEXT)""")
    conn.execute("""CREATE TABLE gpkg_contents (
        table_name TEXT NOT NULL PRIMARY KEY, data_type TEXT NOT NULL DEFAULT 'features',
        identifier TEXT UNIQUE, description TEXT DEFAULT '',
        last_change TEXT NOT NULL DEFAULT '2026-02-10T00:00:00.000Z',
        min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE, srs_id INTEGER)""")
    conn.execute("""CREATE TABLE gpkg_geometry_columns (
        table_name TEXT NOT NULL, column_name TEXT NOT NULL,
        geometry_type_name TEXT NOT NULL, srs_id INTEGER NOT NULL,
        z TINYINT NOT NULL, m TINYINT NOT NULL,
        CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name))""")

    conn.execute("""INSERT INTO gpkg_spatial_ref_sys VALUES
        ('WGS 84', 4326, 'EPSG', 4326,
        'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]', NULL)""")
    conn.execute("INSERT INTO gpkg_spatial_ref_sys VALUES ('Undefined cartesian SRS', -1, 'NONE', -1, 'undefined', NULL)")
    conn.execute("INSERT INTO gpkg_spatial_ref_sys VALUES ('Undefined geographic SRS', 0, 'NONE', 0, 'undefined', NULL)")

    col_defs = ['fid INTEGER PRIMARY KEY AUTOINCREMENT']
    for name, ctype in non_fid_cols:
        col_defs.append(f'"{name}" {ctype}')
    conn.execute(f'CREATE TABLE features ({", ".join(col_defs)})')

    # Envelope
    geom_idx = col_names.index('geom')
    envs = []
    for r in all_rows:
        geom = r[geom_idx]
        if geom and len(geom) > 40:
            envs.append(struct.unpack_from('<dddd', geom, 8))
    minx = min(e[0] for e in envs) if envs else 0
    maxx = max(e[1] for e in envs) if envs else 0
    miny = min(e[2] for e in envs) if envs else 0
    maxy = max(e[3] for e in envs) if envs else 0

    conn.execute("""INSERT INTO gpkg_contents
        (table_name, data_type, identifier, description, min_x, min_y, max_x, max_y, srs_id)
        VALUES ('features', 'features', 'Java_merged', 'Merged East/Central/West Java', ?, ?, ?, ?, 4326)""",
        (minx, miny, maxx, maxy))
    conn.execute("INSERT INTO gpkg_geometry_columns VALUES ('features', 'geom', 'POLYGON', 4326, 0, 0)")

    insert_cols = [f'"{n}"' for n in col_names]
    placeholders = ', '.join(['?'] * len(insert_cols))
    insert_sql = f'INSERT INTO features ({", ".join(insert_cols)}) VALUES ({placeholders})'
    for r in all_rows:
        conn.execute(insert_sql, r)
    conn.commit()

    # Spatial index
    conn.execute('CREATE TABLE rtree_features_geom (id INTEGER PRIMARY KEY, minx REAL, maxx REAL, miny REAL, maxy REAL)')
    for fid, geom in conn.execute('SELECT fid, geom FROM features'):
        if geom and len(geom) > 40:
            env = struct.unpack_from('<dddd', geom, 8)
            conn.execute('INSERT INTO rtree_features_geom VALUES (?,?,?,?,?)',
                         (fid, env[0], env[1], env[2], env[3]))
    conn.commit()

    n = conn.execute('SELECT COUNT(*) FROM features').fetchone()[0]
    conn.close()
    return n


# ─── Main ────────────────────────────────────────────────────────────────────

def run_conversion(base_dir: Path) -> str:
    """Run KML→GPKG conversion. Returns log output as a string."""
    buf = io.StringIO()

    def log(msg: str = "") -> None:
        buf.write(msg + "\n")

    kml_files = sorted(base_dir.glob("*.kml"))

    if not kml_files:
        log("No KML files found in data/kml/")
        log("Place .kml files in the data/kml/ directory and run again.")
        return buf.getvalue()

    log(f"Found {len(kml_files)} KML files\n")

    output_dir = base_dir / "gpkg_output"
    output_dir.mkdir(exist_ok=True)

    # ── Step 1: Convert all KML → GPKG ──
    for kml_path in kml_files:
        name = kml_path.stem
        clean_name = RENAME_MAP.get(name, name)
        gpkg_path = output_dir / f"{clean_name}.gpkg"

        log(f"Processing: {kml_path.name}")
        try:
            features = parse_kml(str(kml_path))
            if not features:
                log(f"  No features found, skipping.\n")
                continue

            n_rows, text_fields = create_gpkg(str(gpkg_path), features, clean_name)
            log(f"  -> {n_rows} features -> {gpkg_path.name}")
            log(f"    Fields: {', '.join(text_fields)} + Longitude, Latitude, Area_Ha\n")

        except Exception as e:
            log(f"  Error: {e}\n")

    # ── Step 2: Post-process the 3 Java files ──
    log("=" * 60)
    log("Post-processing Java province files\n")

    java_gpkgs = []
    for kml_path in kml_files:
        name = kml_path.stem
        if name not in JAVA_STEMS:
            continue
        province = RENAME_MAP[name]
        gpkg_path = output_dir / f"{province}.gpkg"
        if not gpkg_path.exists():
            continue

        log(f"  Post-processing: {province}.gpkg")
        # Redirect prints from postprocess_java_file into our buffer
        old_stdout = sys.stdout
        sys.stdout = buf
        postprocess_java_file(gpkg_path, province, kml_path)
        sys.stdout = old_stdout
        log()

    # ── Step 3: Create merged GPKG (Java files only) ──
    log("=" * 60)
    log("Creating merged GeoPackage: ALL_merged.gpkg (East/Central/West Java only)\n")

    if java_gpkgs:
        merged_path = output_dir / "ALL_merged.gpkg"
        n = create_merged_gpkg(merged_path, java_gpkgs)
        log(f"  -> {n} total features -> ALL_merged.gpkg")

    log("\nDone!")
    log(f"\nOutput directory: {output_dir}")
    return buf.getvalue()


# ─── Team Scripts Plugin ──────────────────────────────────────────────────────


class KmlToGpkgScript:
    name = "KML to GeoPackage"
    description = "Convert uploaded KML files to GeoPackage format with attribute extraction and Java province post-processing."

    async def run(self) -> ScriptResult:
        import zipfile
        try:
            from app.config import KML_DIR
            KML_DIR.mkdir(parents=True, exist_ok=True)
            kml_files = list(KML_DIR.glob("*.kml"))
            if not kml_files:
                return ScriptResult(
                    success=False,
                    output="",
                    error="No KML files found. Upload .kml files first.",
                )
            output = run_conversion(KML_DIR)
            # Zip all .gpkg files for download
            gpkg_dir = KML_DIR / "gpkg_output"
            gpkg_files = sorted(gpkg_dir.glob("*.gpkg")) if gpkg_dir.exists() else []
            zip_path = None
            if gpkg_files:
                zip_path = str(KML_DIR / "geopackages.zip")
                with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                    for gf in gpkg_files:
                        zf.write(gf, gf.name)
                output += f"\n{len(gpkg_files)} GeoPackage file(s) ready for download."
            return ScriptResult(
                success=True,
                output=output,
                download_file=zip_path,
            )
        except Exception as e:
            return ScriptResult(success=False, output="", error=str(e))


script = KmlToGpkgScript()


if __name__ == "__main__":
    from app.config import KML_DIR
    KML_DIR.mkdir(parents=True, exist_ok=True)
    print(run_conversion(KML_DIR))
