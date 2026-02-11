from dataclasses import dataclass
from typing import Protocol


@dataclass
class ScriptResult:
    success: bool
    output: str
    error: str | None = None


class Script(Protocol):
    name: str
    description: str

    async def run(self) -> ScriptResult: ...
