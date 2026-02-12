from pydantic import BaseModel


class ScriptInfo(BaseModel):
    id: str
    name: str
    description: str


class RunResponse(BaseModel):
    success: bool
    output: str
    error: str | None = None
    download_url: str | None = None
