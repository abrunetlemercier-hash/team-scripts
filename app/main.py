from contextlib import asynccontextmanager
from pathlib import Path
from typing import List

from fastapi import FastAPI, Request, UploadFile, File
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from app.config import GPKG_DIR, KML_DIR
from app.models import RunResponse, ScriptInfo
from app.registry import discover_scripts, get_all_scripts, get_script

BASE_DIR = Path(__file__).resolve().parent
INPUT_CONFIG = {
    "kml": {
        "dir": KML_DIR,
        "extension": ".kml",
        "label": "KML Files",
        "empty": "No KML files uploaded yet.",
        "drop_text": "Drop .kml files here or click to browse",
    },
    "gpkg": {
        "dir": GPKG_DIR,
        "extension": ".gpkg",
        "label": "GPKG Files",
        "empty": "No GPKG files uploaded yet.",
        "drop_text": "Drop .gpkg files here or click to browse",
    },
}


@asynccontextmanager
async def lifespan(app: FastAPI):
    KML_DIR.mkdir(parents=True, exist_ok=True)
    GPKG_DIR.mkdir(parents=True, exist_ok=True)
    discover_scripts()
    yield


app = FastAPI(title="Apolownia Tools", lifespan=lifespan)

templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))
app.mount("/static", StaticFiles(directory=str(BASE_DIR / "static")), name="static")


def tool_context(script_id: str, script):
    input_type = getattr(script, "input_type", "kml")
    config = INPUT_CONFIG[input_type]
    files = sorted(f.name for f in config["dir"].glob(f"*{config['extension']}"))
    return {
        "id": script_id,
        "name": script.name,
        "description": script.description,
        "input_type": input_type,
        "input_label": config["label"],
        "input_extension": config["extension"],
        "drop_text": config["drop_text"],
        "empty_text": config["empty"],
        "files": files,
        "has_files": bool(files),
    }


# --- HTML routes ---


@app.get("/", response_class=HTMLResponse)
async def homepage(request: Request):
    scripts = get_all_scripts()
    script_list = []
    for sid, s in scripts.items():
        script_list.append(tool_context(sid, s))
    return templates.TemplateResponse(
        "index.html",
        {"request": request, "scripts": script_list},
    )


@app.get("/tools/{script_id}", response_class=HTMLResponse)
async def tool_page(request: Request, script_id: str):
    script = get_script(script_id)
    if not script:
        return templates.TemplateResponse(
            "result.html",
            {
                "request": request,
                "script_name": script_id,
                "result": None,
                "error": "Tool not found",
                "back_url": "/",
            },
            status_code=404,
        )
    return templates.TemplateResponse(
        "tool.html",
        {"request": request, "tool": tool_context(script_id, script)},
    )


@app.post("/tools/{script_id}/upload")
async def upload_tool_files(script_id: str, files: List[UploadFile] = File(...)):
    script = get_script(script_id)
    if not script:
        return RedirectResponse(url="/", status_code=303)

    input_type = getattr(script, "input_type", "kml")
    config = INPUT_CONFIG[input_type]
    input_dir = config["dir"]
    extension = config["extension"]
    input_dir.mkdir(parents=True, exist_ok=True)

    for f in files:
        if f.filename and f.filename.lower().endswith(extension):
            dest = input_dir / Path(f.filename).name
            content = await f.read()
            dest.write_bytes(content)
    return RedirectResponse(url=f"/tools/{script_id}", status_code=303)


@app.post("/tools/{script_id}/delete/{filename}")
async def delete_tool_file(script_id: str, filename: str):
    script = get_script(script_id)
    if not script:
        return RedirectResponse(url="/", status_code=303)

    input_type = getattr(script, "input_type", "kml")
    config = INPUT_CONFIG[input_type]
    filepath = config["dir"] / Path(filename).name
    if filepath.exists() and filepath.suffix.lower() == config["extension"]:
        filepath.unlink()
    return RedirectResponse(url=f"/tools/{script_id}", status_code=303)


@app.post("/tools/{script_id}/run", response_class=HTMLResponse)
async def run_script_html(request: Request, script_id: str):
    script = get_script(script_id)
    if not script:
        return templates.TemplateResponse(
            "result.html",
            {
                "request": request,
                "script_name": script_id,
                "result": None,
                "error": "Tool not found",
                "back_url": "/",
            },
            status_code=404,
        )
    result = await script.run()
    download_url = None
    if result.download_file:
        # Store the path in a simple in-memory map keyed by filename
        fname = Path(result.download_file).name
        app.state.pending_downloads = getattr(app.state, "pending_downloads", {})
        app.state.pending_downloads[fname] = result.download_file
        download_url = f"/download/{fname}"
    return templates.TemplateResponse(
        "result.html",
        {
            "request": request,
            "script_name": script.name,
            "result": result,
            "download_url": download_url,
            "back_url": f"/tools/{script_id}",
        },
    )


@app.post("/run/{script_id}", response_class=HTMLResponse)
async def run_script_html_legacy(request: Request, script_id: str):
    return await run_script_html(request, script_id)


@app.get("/download/{filename}")
async def download_file(filename: str):
    downloads = getattr(app.state, "pending_downloads", {})
    filepath = downloads.get(filename)
    if not filepath or not Path(filepath).exists():
        return HTMLResponse("File not found", status_code=404)
    return FileResponse(filepath, filename=filename, media_type="application/zip")


@app.post("/refresh")
async def refresh_scripts():
    discover_scripts()
    return RedirectResponse(url="/", status_code=303)


# --- JSON API routes ---


@app.get("/api/scripts")
async def list_scripts_api() -> list[ScriptInfo]:
    scripts = get_all_scripts()
    return [
        ScriptInfo(id=sid, name=s.name, description=s.description)
        for sid, s in scripts.items()
    ]


@app.post("/api/run/{script_id}")
async def run_script_api(script_id: str) -> RunResponse:
    script = get_script(script_id)
    if not script:
        return RunResponse(success=False, output="", error="Script not found")
    result = await script.run()
    return RunResponse(
        success=result.success, output=result.output, error=result.error
    )
