from dataclasses import dataclass
from typing import Protocol


@dataclass
class ScriptResult:
    success: bool
    output: str
    error: str | None = None
    download_file: str | None = None  # absolute path to a file to offer for download


class Script(Protocol):
    name: str
    description: str

    async def run(self) -> ScriptResult: ...
