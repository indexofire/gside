"""sourmash gather-based species identification (species-id plan B)."""

from __future__ import annotations

import csv
import json
import logging
import subprocess
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from gside.analysis._envelope import MethodEnvelope
from gside.analysis._thresholds import ANI_HIGH, ANI_MEDIUM
from gside.analysis.species_name import UNKNOWN, apply_species_name, canonical_species
from gside.config import SPECIES_DB_DIR

# Fallback bands for gather CSVs without the containment-ANI column
# (sourmash >= 4.x always emits it; older exports may not).
_FUW_HIGH = 0.90
_FUW_MEDIUM = 0.70
_MIXTURE_MEMBER = 0.10

logger = logging.getLogger(__name__)


@dataclass
class SourmashIdResult(MethodEnvelope):
    pass


def _float(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "") or 0.0)
    except ValueError:
        return 0.0


def _confidence(ani: float, fuw: float) -> str:
    """Containment ANI when available, hash-fraction bands otherwise."""
    if ani > 0:
        if ani >= ANI_HIGH:
            return "high"
        if ani >= ANI_MEDIUM:
            return "medium"
        return "low"
    if fuw >= _FUW_HIGH:
        return "high"
    if fuw >= _FUW_MEDIUM:
        return "medium"
    return "low"


def _run(cmd: list[str]) -> str:
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if result.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {result.stderr.strip()}")
    return result.stdout


def _gather_rows_legacy(files: Sequence[str | Path], db_zip: Path) -> list[dict[str, Any]]:
    with tempfile.TemporaryDirectory(prefix="gside-sourmash-") as tmp:
        work = Path(tmp)
        sig = work / "query.sig"
        gather_csv = work / "gather.csv"
        _run(
            [
                "sourmash",
                "sketch",
                "dna",
                "-p",
                "k=31,scaled=1000",
                *[str(f) for f in files],
                "-o",
                str(sig),
            ]
        )
        try:
            _run(["sourmash", "gather", str(sig), str(db_zip), "-o", str(gather_csv)])
        except RuntimeError:
            return []
        if not gather_csv.is_file():
            return []
        return list(csv.DictReader(gather_csv.read_text().splitlines()))


def _gather_and_tax(
    contigs: str | Path | Sequence[str | Path], db_dir: Path
) -> list[dict[str, Any]]:
    from gside.engine.gseek_bridge import engine, log_sourmash_shadow, sourmash_gather

    db_zip = db_dir / "gtdb-reps-k31.zip"
    lineages = db_dir / "lineages.csv"
    if not db_zip.exists() or not lineages.exists():
        raise RuntimeError(
            f"sourmash GTDB database not found in {db_dir} (run: gside db setup --tier sourmash)"
        )
    files = [contigs] if isinstance(contigs, (str, Path)) else list(contigs)

    def compare(old: list[dict[str, Any]], new: list[dict[str, Any]]) -> bool:
        log_sourmash_shadow(old, new)
        return True

    rows = engine().run(
        lambda: _gather_rows_legacy(files, db_zip),
        lambda: sourmash_gather(files, db_zip),
        compare=compare,
    )

    lineage_map = _lineage_map(lineages)
    out: list[dict[str, Any]] = []
    for row in rows:
        name = row.get("name", "")
        # gather's `name` is the signature name ("GCF_... Species"); the
        # lineage table is keyed by accession.
        accession = name.split()[0] if name else ""
        species = lineage_map.get(accession) or lineage_map.get(name) or name
        ani = _float(row, "average_containment_ani") or _float(row, "match_containment_ani")
        out.append(
            {
                "species": species,
                "f_unique_weighted": _float(row, "f_unique_weighted"),
                "ani": round(ani * 100.0, 2),
            }
        )
    out.sort(key=lambda hit: hit["f_unique_weighted"], reverse=True)
    return out


def _lineage_map(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    mapping: dict[str, str] = {}
    with path.open() as fh:
        for row in csv.reader(fh):
            if len(row) >= 2:
                mapping[row[0]] = row[-1].split(";")[-1].strip()
    return mapping


def _db_version(db_dir: Path) -> str:
    manifest = db_dir.parent / "manifests" / "D4_sourmash.json"
    if manifest.is_file():
        try:
            checksum = json.loads(manifest.read_text()).get("checksum", "")
            return checksum[:8] if checksum else "RS226"
        except json.JSONDecodeError:
            return "RS226"
    return "RS226"


def identify_by_sourmash(
    contigs: str | Path | Sequence[str | Path], db_dir: str | Path | None = None
) -> SourmashIdResult:
    base = Path(db_dir) if db_dir else SPECIES_DB_DIR / "D4_sourmash"
    partition = _gather_and_tax(contigs, base)
    database = {"name": "D4_sourmash", "version": _db_version(base)}

    best = partition[0] if partition else None
    best_fuw = best["f_unique_weighted"] if best else 0.0
    best_ani = best["ani"] if best else 0.0

    strong_members = [
        p
        for p in partition
        if p["f_unique_weighted"] >= _MIXTURE_MEMBER and (p["ani"] >= ANI_MEDIUM or p["ani"] == 0.0)
    ]
    # Several strains of one species are not a mixture; distinct species are.
    strong_species = {canonical_species(p["species"]) for p in strong_members} - {UNKNOWN}
    flags = ["possible_mixture"] if len(strong_species) >= 2 else []

    confidence = _confidence(best_ani, best_fuw)
    if flags and confidence == "high":
        confidence = "medium"
    species = best["species"] if best is not None and confidence != "low" else "Mixed/Unknown"

    result = {
        "confidence": confidence,
        "f_unique_weighted": best_fuw,
        "ani": best_ani if best_ani > 0 else None,
        "gather_partition": [
            {
                "species": p["species"],
                "f_unique_weighted": p["f_unique_weighted"],
                "ani": p["ani"] if p["ani"] > 0 else None,
            }
            for p in partition
        ],
        "flags": flags,
    }
    apply_species_name(result, species or "Mixed/Unknown")
    return SourmashIdResult(method="sourmash", database=database, result=result)
