"""Database management — download, build, and status for gside reference DBs.

Tiers:
  mini     marker rules + sequences (ships with repo, always available)
  panel    curated reference panel (skani sketch, ~2.7GB)
  mash     RefSeq MinHash sketch (~331MB)
  all      panel + mash
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from .config import DATA_DIR, SPECIES_DB_DIR, which

_TIERS: dict[str, dict[str, Any]] = {
    "mini": {
        "description": "Marker rules + sequences (bundled, always ready)",
        "size": "~1MB",
        "items": ["markers_v2"],
    },
    "panel": {
        "description": "Curated reference panel for ANI (skani)",
        "size": "~2.7GB",
        "items": ["refseq_panel"],
    },
    "mash": {
        "description": "RefSeq MinHash sketch for distance-based ID",
        "size": "~331MB",
        "items": ["mash_refseq"],
    },
    "all": {
        "description": "panel + mash (full multi-layer capability)",
        "size": "~3GB",
        "items": ["refseq_panel", "mash_refseq"],
    },
}


def db_status() -> dict[str, dict[str, Any]]:
    markers_ok = (DATA_DIR / "reference/species/markers_v2.fasta").exists()
    panel_ok = _check_panel()
    mash_ok = _check_mash()

    return {
        "markers_v2": {
            "ready": markers_ok,
            "path": str(DATA_DIR / "reference/species"),
            "tier": "mini",
        },
        "refseq_panel": {
            "ready": panel_ok,
            "path": str(SPECIES_DB_DIR / "refseq_panel"),
            "tier": "panel",
        },
        "mash_refseq": {
            "ready": mash_ok,
            "path": str(SPECIES_DB_DIR / "mash_refseq"),
            "tier": "mash",
        },
    }


def _check_panel() -> bool:
    p = SPECIES_DB_DIR / "refseq_panel"
    return (p / "panel.sketch" / "sketches.db").exists() or (
        p / "panel.sketch"
    ).is_file()


def _check_mash() -> bool:
    m = SPECIES_DB_DIR / "mash_refseq"
    return (m / "mash.msh").exists() or (m / "payload.bin").exists()


def db_setup(
    tier: str = "panel",
    source: str = "",
) -> dict[str, str]:
    if tier not in _TIERS:
        raise ValueError(f"unknown tier {tier!r}; expected one of {sorted(_TIERS)}")

    results: dict[str, str] = {}

    if tier in ("panel", "all"):
        results["refseq_panel"] = _setup_panel(source)
    if tier in ("mash", "all"):
        results["mash_refseq"] = _setup_mash(source)

    if not results:
        results["info"] = "mini tier ships with repo — nothing to download"

    return results


def _setup_panel(source: str) -> str:
    if _check_panel():
        return "already ready"

    if source:
        src = Path(source)
        if src.name != "refseq_panel":
            src = src / "refseq_panel"
        if src.exists():
            dst = SPECIES_DB_DIR / "refseq_panel"
            dst.mkdir(parents=True, exist_ok=True)
            shutil.copytree(src, dst, dirs_exist_ok=True)
            return f"copied from {src}"

    skani = which("skani")
    if not skani:
        return "ERROR: skani not found (pixi add skani)"

    genomes_dir = SPECIES_DB_DIR / "refseq_panel" / "genomes"
    if not genomes_dir.exists() or not list(genomes_dir.glob("*.fna")):
        return "ERROR: no genomes found — run bacmap's curate_species_panel.py first or provide --source"

    sketch = SPECIES_DB_DIR / "refseq_panel" / "panel.sketch"
    result = subprocess.run(
        [skani, "sketch", *[str(f) for f in sorted(genomes_dir.glob("*.fna"))], "-o", str(sketch)],
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if result.returncode != 0:
        return f"ERROR: skani sketch failed: {result.stderr[:200]}"
    return f"built panel.sketch from {len(list(genomes_dir.glob('*.fna')))} genomes"


def _setup_mash(source: str) -> str:
    if _check_mash():
        return "already ready"

    if source:
        src = Path(source)
        if src.name != "mash_refseq":
            src = src / "mash_refseq"
        if src.exists():
            dst = SPECIES_DB_DIR / "mash_refseq"
            dst.mkdir(parents=True, exist_ok=True)
            shutil.copytree(src, dst, dirs_exist_ok=True)
            return f"copied from {src}"

    mash = which("mash")
    if not mash:
        return "ERROR: mash not found (pixi add mash)"

    return "ERROR: mash sketch requires RefSeq FASTAs — provide --source pointing to bacmap's data/db/mash_refseq"


def run_db_command(args: list[str]) -> int:
    if not args:
        _print_status()
        return 0

    sub = args[0]
    if sub == "status":
        _print_status()
        return 0
    elif sub == "setup":
        tier = args[1] if len(args) > 1 else "panel"
        source = ""
        if "--source" in args:
            i = args.index("--source")
            source = args[i + 1] if i + 1 < len(args) else ""
        results = db_setup(tier, source)
        for name, msg in results.items():
            print(f"  {name}: {msg}")
        has_error = any("ERROR" in v for v in results.values())
        return 1 if has_error else 0
    elif sub == "list":
        for name, info in _TIERS.items():
            print(f"  {name:8} {info['size']:8} {info['description']}")
        return 0
    else:
        print(f"unknown subcommand: {sub} (status|setup|list)")
        return 1


def _print_status() -> None:
    status = db_status()
    print("gside database status")
    print("─" * 50)
    for name, info in status.items():
        icon = "✅" if info["ready"] else "❌"
        print(f"  {icon} {name:16} tier={info['tier']:6} path={info['path']}")
    print()
    print("  Run 'gside db setup --help' to provision")
