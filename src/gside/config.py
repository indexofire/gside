"""Runtime configuration: data dirs and pixi-aware binary resolution."""

from __future__ import annotations

import os
import shutil
from collections.abc import Mapping
from enum import StrEnum
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_ROOT.parent.parent


def _default_data_dir(package_root: Path, project_root: Path) -> Path:
    repo_data = project_root / "data"
    return repo_data if repo_data.is_dir() else package_root / "data"


DATA_DIR = Path(os.environ.get("GSIDE_DATA_DIR", _default_data_dir(PACKAGE_ROOT, PROJECT_ROOT)))


class DbRootSource(StrEnum):
    """Where the reference-database root came from (`db status` reports it)."""

    ENV = "env"
    REPO = "repo"
    XDG = "xdg"


ENV_DB_DIR = "GSIDE_DB_DIR"
ENV_DATA_DIR = "GSIDE_DATA_DIR"
ENV_XDG_DATA_HOME = "XDG_DATA_HOME"
#: XDG base-directory default for the database root (gspec datadir convention).
XDG_DB_ROOT = Path("~/.local/share/gside/db")

#: Wheel-bundled database seed (force-included from ``data/db`` at build
#: time). Source checkouts keep their databases in the repo, so this exists
#: only inside installed packages; ``gside.db.ensure_marker_seed`` copies
#: from it into the user's root on first use.
BUNDLED_DB_ROOT = PACKAGE_ROOT / "data" / "db"


def _root_candidates(
    environ: Mapping[str, str],
    repo_data_dir: Path | None,
) -> list[tuple[DbRootSource, Path]]:
    """Active (source, path) candidates in precedence order (never empty)."""
    candidates: list[tuple[DbRootSource, Path]] = []
    raw_db = environ.get(ENV_DB_DIR, "").strip()
    raw_data = environ.get(ENV_DATA_DIR, "").strip()
    if raw_db:
        candidates.append((DbRootSource.ENV, Path(raw_db).expanduser()))
    elif raw_data:
        candidates.append((DbRootSource.ENV, Path(raw_data).expanduser() / "db"))
    repo = PROJECT_ROOT / "data" if repo_data_dir is None else repo_data_dir
    if repo.is_dir():
        candidates.append((DbRootSource.REPO, repo / "db"))
    xdg = environ.get(ENV_XDG_DATA_HOME, "").strip()
    xdg_path = Path(xdg).expanduser() if xdg else None
    if xdg_path is not None and xdg_path.is_absolute():
        candidates.append((DbRootSource.XDG, xdg_path / "gside" / "db"))
    else:
        candidates.append((DbRootSource.XDG, XDG_DB_ROOT.expanduser()))
    return candidates


def resolve_db_root(
    environ: Mapping[str, str] | None = None,
    *,
    repo_data_dir: Path | None = None,
) -> tuple[Path, DbRootSource]:
    """Resolve the reference-database root and where the choice came from.

    Precedence: ``$GSIDE_DB_DIR`` > ``$GSIDE_DATA_DIR/db`` > the source
    checkout's ``<repo>/data/db`` > ``$XDG_DATA_HOME/gside/db`` >
    ``~/.local/share/gside/db``. A relative or empty ``$XDG_DATA_HOME`` is
    ignored per the XDG Base Directory specification. The returned path is
    not required to exist; ``gside.db.ensure_db_root`` creates it on write.
    """
    env = os.environ if environ is None else environ
    source, path = _root_candidates(env, repo_data_dir)[0]
    return path, source


def db_root_source(
    path: Path | None = None,
    environ: Mapping[str, str] | None = None,
    *,
    repo_data_dir: Path | None = None,
) -> str:
    """Classify *path* (default: the effective root) as env/repo/xdg/custom.

    ``custom`` covers roots that match no active candidate — a test-injected
    directory or code that reassigns ``SPECIES_DB_DIR`` after import.
    """
    env = os.environ if environ is None else environ
    target = (SPECIES_DB_DIR if path is None else Path(path)).expanduser().resolve()
    for source, candidate in _root_candidates(env, repo_data_dir):
        if candidate.expanduser().resolve() == target:
            return source.value
    return "custom"


SPECIES_DB_DIR, DB_ROOT_SOURCE = resolve_db_root()


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
