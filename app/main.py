from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from app.models import RunResponse, ScriptInfo
from app.registry import discover_scripts, get_all_scripts, get_script

BASE_DIR = Path(__file__).resolve().parent


@asynccontextmanager
async def lifespan(app: FastAPI):
    discover_scripts()
    yield


app = FastAPI(title="Team Scripts", lifespan=lifespan)

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
    return templates.TemplateResponse(
        "index.html", {"request": request, "scripts": script_list}
    )


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
    return templates.TemplateResponse(
        "result.html",
        {"request": request, "script_name": script.name, "result": result},
    )


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
