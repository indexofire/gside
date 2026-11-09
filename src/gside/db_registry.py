"""Named per-mode databases: resolution, enumeration, and validation.

A mode can draw on several databases. Besides the bundled/setup default tier
(``D1_marker``/``D2_ani``/``D3_mash``/``D4_sourmash``), named databases live
under the mode's family directory: ``<root>/<family>/<name>`` — e.g.
``data/db/marker/food`` or ``data/db/panel/lab2026``. A ``-d/--db`` value is
resolved in this order:

1. an existing directory path → used directly (layout must validate);
2. the reserved name ``default`` → the mode's tier directory under the root;
3. a name under ``<root>/<family>/``.

``MODE:NAME`` qualifies a selection when several database-capable modes run
in one invocation (``-m marker+panel``). Modes without a family (``gtdbtk``,
meta modes) reject selections with an actionable error.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from gside import config

FAMILY_BY_MODE: dict[str, str] = {
    "marker": "marker",
    "panel": "panel",
    "mash_refseq": "mash_refseq",
    "sourmash": "sourmash",
}

_DEFAULT_DIR_BY_MODE: dict[str, str] = {
    "marker": "D1_marker",
    "panel": "D2_ani",
    "mash_refseq": "D3_mash",
    "sourmash": "D4_sourmash",
}

DEFAULT_NAME = "default"

__all__ = [
    "DbChoice",
    "FAMILY_BY_MODE",
    "DEFAULT_NAME",
    "enumerate_databases",
    "layout_valid",
    "parse_db_value",
    "resolve_db",
    "resolve_selections",
]


@dataclass(frozen=True)
class DbChoice:
    """One resolved database for one mode."""

    mode: str
    name: str
    path: Path
    ready: bool = True


def _effective_root(root: str | Path | None) -> Path:
    return Path(root) if root else config.SPECIES_DB_DIR


def _layout_requirements(mode: str) -> tuple[tuple[str, ...], ...]:
    """Alternative file groups that make a directory valid for *mode*.

    Marker databases keep the stricter ``MarkerDatabase`` contract (checked
    directly in :func:`layout_valid`); the groups below mirror
    ``db._COMPONENTS`` readiness probes for the sketch families.
    """
    if mode == "panel":
        return (("panel.sketch/sketches.db",), ("panel.sketch",))
    if mode == "mash_refseq":
        return (("mash.msh",),)
    if mode == "sourmash":
        return (("gtdb-reps-k31.zip", "lineages.csv"),)
    raise ValueError(f"{mode} does not support named databases")


def layout_valid(mode: str, path: Path) -> bool:
    """True when *path* holds a usable database of *mode*'s family."""
    if not path.is_dir():
        return False
    if mode == "marker":
        from gside.marker_db import MarkerDatabase

        try:
            MarkerDatabase.from_path(path)
        except ValueError:
            return False
        return True
    groups = _layout_requirements(mode)
    return any(all((path / rel).exists() for rel in group) for group in groups)


def parse_db_value(value: str) -> tuple[str | None, str]:
    """``MODE:NAME`` → (mode, name); a bare name or path → (None, value)."""
    raw = value.strip()
    if not raw or raw.endswith(":") or raw.startswith(":"):
        raise ValueError(f"empty database name in {value!r}")
    if ":" in raw:
        mode, _, name = raw.partition(":")
        return mode.strip(), name.strip()
    return None, raw


def _available_names(mode: str, root: Path) -> list[str]:
    family = root / FAMILY_BY_MODE[mode]
    if not family.is_dir():
        return []
    return sorted(
        child.name
        for child in family.iterdir()
        if child.is_dir() and layout_valid(mode, child) and not child.name.startswith(".")
    )


def resolve_db(mode: str, value: str, root: str | Path | None = None) -> DbChoice:
    """Resolve one ``-d`` value for *mode* into a concrete database directory."""
    if mode not in FAMILY_BY_MODE:
        raise ValueError(
            f"{mode} does not support named databases; "
            "pick one of: marker, panel, mash_refseq, sourmash"
        )
    base = _effective_root(root)
    raw = value.strip()

    candidate = Path(raw).expanduser()
    if candidate.is_dir():
        if not layout_valid(mode, candidate):
            raise ValueError(f"{candidate} is not a valid {mode} database")
        return DbChoice(mode=mode, name=candidate.name, path=candidate.resolve())

    if raw == DEFAULT_NAME:
        tier = base / _DEFAULT_DIR_BY_MODE[mode]
        return DbChoice(mode=mode, name=DEFAULT_NAME, path=tier, ready=layout_valid(mode, tier))

    named = base / FAMILY_BY_MODE[mode] / raw
    if named.is_dir() and not layout_valid(mode, named):
        raise ValueError(f"{named} is not a valid {mode} database")
    if not layout_valid(mode, named):
        available = _available_names(mode, base)
        hint = f"; available: {', '.join(available)}" if available else ""
        raise ValueError(f"no {mode} database named {raw!r} under {base}{hint}")
    return DbChoice(mode=mode, name=raw, path=named)


def resolve_selections(
    values: list[str] | tuple[str, ...],
    modes: list[str] | tuple[str, ...],
    root: str | Path | None = None,
) -> dict[str, list[DbChoice]]:
    """Resolve every ``-d`` value against the invocation's mode list.

    Unqualified names need exactly one database-capable mode in the
    invocation; qualified names must name one of the invocation's modes.
    Duplicates collapse in first-seen order.
    """
    active_family_modes = [m for m in modes if m in FAMILY_BY_MODE]
    resolved: dict[str, list[DbChoice]] = {}
    seen: dict[str, list[Path]] = {}
    for value in values:
        mode, name = parse_db_value(value)
        if mode is None and len(active_family_modes) == 1:
            mode = active_family_modes[0]
        if mode is None:
            if not active_family_modes:
                raise ValueError(
                    "no mode in this invocation supports named databases "
                    "(supported: marker, panel, mash_refseq, sourmash)"
                )
            raise ValueError(
                "several database-capable modes in this invocation "
                f"({', '.join(active_family_modes)}); qualify the database: "
                "e.g. --db marker:food --db panel:lab2026"
            )
        if mode not in modes:
            raise ValueError(
                f"--db {value!r} names mode {mode!r}, which is not part of this invocation "
                f"({', '.join(modes)})"
            )
        choice = resolve_db(mode, name, root=root)
        known = seen.setdefault(mode, [])
        if choice.path not in known:
            known.append(choice.path)
            resolved.setdefault(mode, []).append(choice)
    return resolved


def enumerate_databases(root: str | Path | None = None) -> dict[str, list[DbChoice]]:
    """Default tier + every valid named database, per database-capable mode."""
    base = _effective_root(root)
    databases: dict[str, list[DbChoice]] = {}
    for mode in FAMILY_BY_MODE:
        tier = base / _DEFAULT_DIR_BY_MODE[mode]
        entries = [
            DbChoice(mode=mode, name=DEFAULT_NAME, path=tier, ready=layout_valid(mode, tier))
        ]
        entries.extend(
            DbChoice(mode=mode, name=name, path=base / FAMILY_BY_MODE[mode] / name)
            for name in _available_names(mode, base)
        )
        databases[mode] = entries
    return databases
