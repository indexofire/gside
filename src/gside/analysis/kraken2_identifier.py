"""Kraken2 + Bracken metagenomic profiling (gside meta, read/contig based)."""

from __future__ import annotations

import tempfile
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from gside.config import SPECIES_DB_DIR, which
from gside.engine._env import require_bin, run_checked

_K2_HINT = "kraken2 not found — pixi add kraken2 bracken"
_BRACKEN_HINT = "bracken not found — pixi add kraken2 bracken"
_K2_INDEXES_URL = "https://benlangmead.github.io/aws-indexes/k2"

_HIGH_FRACTION = 0.10
_MEDIUM_FRACTION = 0.01
_TOP_N = 20
_BRACKEN_READ_LEN = "150"
_BRACKEN_THRESHOLD = "10"


@dataclass
class Kraken2Result:
    method: str = "kraken2_bracken"
    database: dict[str, str] = field(default_factory=dict)
    result: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "analysis_type": "metagenomic_profiling",
            "method": self.method,
            "database": self.database,
            "result": self.result,
        }


def _fraction_confidence(fraction: float) -> str:
    if fraction >= _HIGH_FRACTION:
        return "high"
    if fraction >= _MEDIUM_FRACTION:
        return "medium"
    return "low"


def _parse_bracken(path: Path) -> list[dict[str, Any]]:
    """Species rows of a bracken report, sorted by fraction descending.

    bracken emits kraken-style rank codes ("S" for species); the spelled-out
    "species" is accepted too so hand-written reports keep parsing.
    """
    species: list[dict[str, Any]] = []
    lines = path.read_text().splitlines()
    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) < 7:
            continue
        name, taxid, lvl, _assigned, _added, new_est, fraction = cols[:7]
        if lvl not in ("S", "species"):
            continue
        try:
            taxid = int(taxid)
            reads = int(new_est)
            frac = float(fraction)
        except ValueError:
            continue
        species.append(
            {
                "name": name,
                "taxid": taxid,
                "reads": reads,
                "fraction": frac,
                "confidence": _fraction_confidence(frac),
            }
        )
    species.sort(key=lambda s: s["fraction"], reverse=True)
    return species


def _parse_kreport(path: Path) -> tuple[int, int]:
    """(total_reads, classified_reads) from the root (R) and unclassified (U) rows."""
    classified = unclassified = 0
    for line in path.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 5:
            continue
        code = cols[3]
        if code == "R":
            classified = int(cols[1])
        elif code == "U":
            unclassified = int(cols[1])
    return classified + unclassified, classified


def identify_by_kraken2(
    files: Sequence[str],
    db_dir: str | Path | None = None,
    threads: int = 4,
) -> Kraken2Result:
    """Profile one sample with kraken2 + bracken (1 file = single-end/contigs, 2 = paired)."""
    k2_bin = require_bin("kraken2", hint=_K2_HINT, resolver=which)
    bracken_bin = require_bin("bracken", hint=_BRACKEN_HINT, resolver=which)

    db = Path(db_dir) if db_dir else SPECIES_DB_DIR / "D5_kraken2"
    if not any(db.glob("*.k2d")):
        raise RuntimeError(
            f"Kraken2 database not found in {db} (download a prebuilt index: {_K2_INDEXES_URL})"
        )

    inputs = [str(f) for f in files]
    with tempfile.TemporaryDirectory(prefix="gside-kraken2-") as tmp:
        work = Path(tmp)
        kreport = work / "kreport.txt"
        kraken_cmd = [
            k2_bin,
            "--db",
            str(db),
            "--threads",
            str(threads),
            "--output",
            str(work / "kraken.out"),
            "--report",
            str(kreport),
        ]
        if len(inputs) == 2:
            kraken_cmd.append("--paired")
        kraken_cmd += inputs
        run_checked(kraken_cmd, timeout=3600, name="kraken2 failed")

        bracken_out = work / "bracken.txt"
        run_checked(
            [
                bracken_bin,
                "-d",
                str(db),
                "-i",
                str(kreport),
                "-o",
                str(bracken_out),
                "-r",
                _BRACKEN_READ_LEN,
                "-t",
                _BRACKEN_THRESHOLD,
            ],
            timeout=300,
            name="bracken failed",
        )
        species = _parse_bracken(bracken_out)
        total_reads, classified_reads = _parse_kreport(kreport)

    return Kraken2Result(
        database={"name": "D5_kraken2", "path": str(db)},
        result={
            "top_species": species[:_TOP_N],
            "total_reads": total_reads,
            "classified_reads": classified_reads,
            "n_species": len(species),
        },
    )
