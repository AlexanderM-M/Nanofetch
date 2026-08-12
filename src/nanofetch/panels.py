from importlib import resources
from typing import Dict, List

from .errors import NanoFetchError


PANELS_PACKAGE = "nanofetch.data.panels"


def available_panels() -> List[str]:
    return sorted(path.name[:-4] for path in resources.files(PANELS_PACKAGE).iterdir()
                  if path.name.endswith(".txt"))


def load_panel(name: str) -> List[str]:
    normalized = name.lower()
    if normalized not in available_panels():
        choices = ", ".join(available_panels()) or "none"
        raise NanoFetchError(f"Unknown panel {name!r}. Available panels: {choices}.")
    resource = resources.files(PANELS_PACKAGE).joinpath(f"{normalized}.txt")
    with resource.open("r", encoding="utf-8") as handle:
        return [line.strip() for line in handle
                if line.strip() and not line.lstrip().startswith("#")]


def panel_descriptions() -> Dict[str, str]:
    return {"cns": "Convenience CNS tumor gene set (not a clinical assay)"}

