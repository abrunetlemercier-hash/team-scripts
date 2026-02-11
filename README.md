# Team Scripts

Web UI for running team automation scripts. Built with FastAPI.

## Quick Start

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload
```

Open http://localhost:8000

## Adding a Script

Create a new file in `scripts/`, e.g. `scripts/my_script.py`:

```python
from scripts._base import ScriptResult

class MyScript:
    name = "My Script"
    description = "Does something useful."

    async def run(self) -> ScriptResult:
        # Your logic here
        return ScriptResult(success=True, output="Done!")

script = MyScript()
```

The script appears automatically on the homepage. No other files to edit.

## API

- `GET /api/scripts` - List all scripts (JSON)
- `POST /api/run/{script_id}` - Run a script (JSON)

## Deploy to Render

1. Push this repo to GitHub
2. Create a new **Web Service** on [render.com](https://render.com)
3. Connect the GitHub repo (Render detects the Dockerfile)
4. Select the free instance type
5. Deploy — share the `.onrender.com` URL with your team
