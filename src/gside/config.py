"""Runtime configuration: data dirs and pixi-aware binary resolution."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_ROOT.parent.parent


def _default_data_dir(package_root: Path, project_root: Path) -> Path:
    repo_data = project_root / "data"
    return repo_data if repo_data.is_dir() else package_root / "data"


DATA_DIR = Path(os.environ.get("GSIDE_DATA_DIR", _default_data_dir(PACKAGE_ROOT, PROJECT_ROOT)))
SPECIES_DB_DIR = Path(os.environ.get("GSIDE_DB_DIR", DATA_DIR / "db"))


def _env_path(*names: str) -> Path | None:
    for name in names:
        value = os.environ.get(name, "")
        if value:
            return Path(value)
    return None


CHECKM2_DB = _env_path("CHECKM2DB")
GTDB_DB = _env_path("GTDBTK_DATA_PATH", "GTDBDB")


# Search order: GSIDE_PIXI_BIN, the repo's own pixi env, the user-level pixi
# env, then PATH. Machine-specific tool dirs belong in GSIDE_PIXI_BIN.
_KNOWN_PIXI_ROOTS = (
    PROJECT_ROOT / ".pixi/envs/default/bin",
    Path.home() / ".pixi/envs/default/bin",
)


def pixi_path() -> str:
    candidates = [os.environ.get("GSIDE_PIXI_BIN", "")]
    candidates += [str(p) for p in _KNOWN_PIXI_ROOTS if p.exists()]
    prefix = ":".join(c for c in candidates if c)
    return f"{prefix}:{os.environ.get('PATH', '')}"


def which(tool: str) -> str | None:
    return shutil.which(tool, path=pixi_path())
