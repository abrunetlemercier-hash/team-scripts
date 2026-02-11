import asyncio
import platform
from datetime import datetime, timezone

from scripts._base import ScriptResult


class SystemInfoScript:
    name = "System Info"
    description = "Collects basic system information (hostname, OS, Python version)."

    async def run(self) -> ScriptResult:
        try:
            lines = [
                f"Hostname:       {platform.node()}",
                f"OS:             {platform.system()} {platform.release()}",
                f"Architecture:   {platform.machine()}",
                f"Python:         {platform.python_version()}",
                f"Timestamp:      {datetime.now(timezone.utc).isoformat()}",
            ]

            proc = await asyncio.create_subprocess_exec(
                "uptime",
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, _ = await proc.communicate()
            if stdout:
                lines.append(f"Uptime:         {stdout.decode().strip()}")

            return ScriptResult(success=True, output="\n".join(lines))
        except Exception as e:
            return ScriptResult(success=False, output="", error=str(e))


script = SystemInfoScript()
