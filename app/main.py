from contextlib import asynccontextmanager
from pathlib import Path
from typing import List

from fastapi import FastAPI, Request, UploadFile, File
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from app.config import KML_DIR, PROJECT_ROOT
from app.models import RunResponse, ScriptInfo
from app.registry import discover_scripts, get_all_scripts, get_script

BASE_DIR = Path(__file__).resolve().parent


@asynccontextmanager
async def lifespan(app: FastAPI):
    KML_DIR.mkdir(parents=True, exist_ok=True)
    discover_scripts()
    yield


app = FastAPI(title="Apolownia Tools", lifespan=lifespan)

templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))
app.mount("/static", StaticFiles(directory=str(BASE_DIR / "static")), name="static")


# --- HTML routes ---


@app.get("/", response_class=HTMLResponse)
async def homepage(request: Request):
    scripts = get_all_scripts()
    script_list = [
        {"id": sid, "name": s.name, "description": s.description}
        for sid, s in scripts.items()
    ]
    kml_files = sorted(f.name for f in KML_DIR.glob("*.kml"))
    return templates.TemplateResponse(
        "index.html",
        {"request": request, "scripts": script_list, "kml_files": kml_files},
    )


@app.post("/upload")
async def upload_kml(files: List[UploadFile] = File(...)):
    KML_DIR.mkdir(parents=True, exist_ok=True)
    for f in files:
        if f.filename and f.filename.lower().endswith(".kml"):
            dest = KML_DIR / f.filename
            content = await f.read()
            dest.write_bytes(content)
    return RedirectResponse(url="/", status_code=303)


@app.post("/delete-kml/{filename}")
async def delete_kml(filename: str):
    filepath = KML_DIR / filename
    if filepath.exists() and filepath.suffix.lower() == ".kml":
        filepath.unlink()
    return RedirectResponse(url="/", status_code=303)


@app.post("/run/{script_id}", response_class=HTMLResponse)
async def run_script_html(request: Request, script_id: str):
    script = get_script(script_id)
    if not script:
        return templates.TemplateResponse(
            "result.html",
            {
                "request": request,
                "script_name": script_id,
                "result": None,
                "error": "Script not found",
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
        },
    )


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
