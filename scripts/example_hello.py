from scripts._base import ScriptResult


class HelloScript:
    name = "Hello World"
    description = "A simple test script that prints a greeting."

    async def run(self) -> ScriptResult:
        return ScriptResult(
            success=True,
            output="Hello from Team Scripts! Everything is working.",
        )


script = HelloScript()
