"""GTDB-Tk phylogenetic species identification (layer 3, most authoritative)."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from gside.config import GTDB_DB, which
from gside.engine._env import require_bin, run_checked

_UNAVAILABLE = (
    "GTDB-Tk not available — install gtdbtk and set GTDBTK_DATA_PATH"
    " (see https://ecogenomics.github.io/GTDBTk/)"
)


@dataclass
class GtdbkResult:
    method: str = "gtdbtk"
    database: dict[str, str] = field(default_factory=dict)
    result: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "analysis_type": "species_identification",
            "method": self.method,
            "database": self.database,
            "result": self.result,
        }


def _species_from_classification(classification: str) -> str:
    """Return the species name from the last ``s__`` rank of a GTDB-Tk lineage."""
    for rank in reversed(classification.split(";")):
        name = rank.strip()
        if name.startswith("s__"):
            species = name[3:].strip()
            return species or "Unknown"
    return "Unknown"


def _read_classification(out_dir: Path) -> str:
    """Parse ``gtdbtk.bac120.summary.tsv`` (fallback ``ar122``) for the classification."""
    for summary_name in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar122.summary.tsv"):
        summary = out_dir / summary_name
        if not summary.is_file():
            continue
        lines = summary.read_text().strip().split("\n")
        if len(lines) < 2:
            continue
        header = lines[0].split("\t")
        data = lines[1].split("\t")
        row = dict(zip(header, data, strict=False))
        return row.get("classification", "")
    return ""


def identify_by_gtdbtk(contigs: str | Path, threads: int = 4) -> GtdbkResult:
    """Classify one assembly with GTDB-Tk ``classify_wf`` (minutes per genome)."""
    if GTDB_DB is None:
        raise RuntimeError(_UNAVAILABLE)
    binary = require_bin("gtdbtk", hint=_UNAVAILABLE, resolver=which)
    with tempfile.TemporaryDirectory(prefix="gside-gtdbtk-") as tmp:
        genome_dir = Path(tmp) / "genomes"
        genome_dir.mkdir()
        (genome_dir / "query.fna").write_text(Path(contigs).read_text())
        out_dir = Path(tmp) / "output"
        run_checked(
            [
                binary,
                "classify_wf",
                "--genome_dir",
                str(genome_dir),
                "--out_dir",
                str(out_dir),
                "--database",
                str(GTDB_DB),
                "--cpus",
                str(threads),
            ],
            timeout=3600,
            name="gtdbtk",
        )
        classification = _read_classification(out_dir)
    species = _species_from_classification(classification)
    result = {
        "species": species,
        "confidence": "high" if species != "Unknown" else "low",
        "classification": classification,
    }
    return GtdbkResult(
        method="gtdbtk",
        database={"name": "GTDB-Tk", "data_path": str(GTDB_DB)},
        result=result,
    )
