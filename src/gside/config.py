"""Runtime configuration: data dirs and pixi-aware binary resolution."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_ROOT.parent.parent

DATA_DIR = Path(os.environ.get("GSIDE_DATA_DIR", PROJECT_ROOT / "data"))
SPECIES_DB_DIR = Path(os.environ.get("GSIDE_DB_DIR", DATA_DIR / "db"))


def _env_path(*names: str) -> Path | None:
    for name in names:
        value = os.environ.get(name, "")
        if value:
            return Path(value)
    return None


CHECKM2_DB = _env_path("CHECKM2DB")
GTDB_DB = _env_path("GTDBTK_DATA_PATH", "GTDBDB")


_KNOWN_PIXI_ROOTS = (
    PROJECT_ROOT / ".pixi/envs/default/bin",
    Path.home() / ".pixi/envs/default/bin",
    Path("/home/mark/repos/github/hermes-bacmap/.pixi/envs/default/bin"),
)


def pixi_path() -> str:
    candidates = [os.environ.get("GSIDE_PIXI_BIN", "")]
    candidates += [str(p) for p in _KNOWN_PIXI_ROOTS if p.exists()]
    prefix = ":".join(c for c in candidates if c)
    return f"{prefix}:{os.environ.get('PATH', '')}"


def which(tool: str) -> str | None:
    return shutil.which(tool, path=pixi_path())
