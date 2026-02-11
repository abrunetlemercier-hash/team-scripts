import importlib
import pkgutil
from pathlib import Path

from scripts._base import Script

_registry: dict[str, Script] = {}


def discover_scripts() -> None:
    """Scan the scripts/ package and register every module that exposes a `script` attribute."""
    _registry.clear()
    scripts_path = Path(__file__).resolve().parent.parent / "scripts"

    for module_info in pkgutil.iter_modules([str(scripts_path)]):
        if module_info.name.startswith("_"):
            continue
        module = importlib.import_module(f"scripts.{module_info.name}")
        if hasattr(module, "script"):
            _registry[module_info.name] = module.script


def get_all_scripts() -> dict[str, Script]:
    return dict(_registry)


def get_script(script_id: str) -> Script | None:
    return _registry.get(script_id)
