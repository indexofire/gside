"""Sylph ANI-based metagenomic profiling (gside meta, FASTQ reads only)."""

from __future__ import annotations

import tempfile
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from gside.config import SPECIES_DB_DIR, which
from gside.engine._env import require_bin, run_checked

_SYLPH_HINT = "sylph not found — pixi add sylph"
_SYLPH_DOCS_URL = "https://sylph-docs.github.io/"
FASTA_INPUT_ERROR = "sylph requires FASTQ reads, not contigs — use --mode kraken2 for FASTA input"

_FASTA_SUFFIXES = (".fna", ".fa", ".fasta")

# Same thresholds as the panel ANI layer (cli.py:_arbitrate high/medium split).
_HIGH_ANI = 95.0
_MEDIUM_ANI = 93.0


@dataclass
class SylphResult:
    method: str = "sylph_ani"
    database: dict[str, str] = field(default_factory=dict)
    result: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "analysis_type": "metagenomic_profiling",
            "method": self.method,
            "database": self.database,
            "result": self.result,
        }


def _ani_confidence(ani: float) -> str:
    if ani >= _HIGH_ANI:
        return "high"
    if ani >= _MEDIUM_ANI:
        return "medium"
    return "low"


def _looks_like_fasta(files: Sequence[str]) -> bool:
    return any(Path(f).suffix in _FASTA_SUFFIXES for f in files)


def _species_from_genome_file(genome_file: str) -> str:
    """`Escherichia_coli.fna` → `Escherichia coli` (filename, extension stripped)."""
    return Path(genome_file).stem.replace("_", " ")


def _find_sketch(work: Path) -> Path:
    """Locate the read sketch sylph wrote (v0 emits .syldata, v1 emits .sylsp
    inside the ``-d`` sketches directory)."""
    candidates = sorted(work.rglob("*.sylsp")) + sorted(work.rglob("*.syldata"))
    if not candidates:
        raise RuntimeError("sylph sketch produced no .syldata/.sylsp output")
    return candidates[0]


def _parse_profile(path: Path) -> list[dict[str, Any]]:
    """sylph profile TSV → species rows sorted by taxonomic abundance descending."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    species: list[dict[str, Any]] = []
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t"), strict=False))
        try:
            ani = float(row.get("Adjusted_ANI", 0.0))
            abundance = float(row.get("Taxonomic_abundance", 0.0))
            coverage = float(row.get("Eff_cov", 0.0))
        except ValueError:
            continue
        species.append(
            {
                "name": _species_from_genome_file(row.get("Genome_file", "")),
                "ani": ani,
                "abundance": abundance,
                "coverage": coverage,
                "confidence": _ani_confidence(ani),
            }
        )
    species.sort(key=lambda s: s["abundance"], reverse=True)
    return species


def identify_by_sylph(
    files: Sequence[str],
    db_dir: str | Path | None = None,
    threads: int = 4,
) -> SylphResult:
    """Profile one FASTQ sample with sylph sketch + profile (1 file = single, 2 = paired)."""
    sylph_bin = require_bin("sylph", hint=_SYLPH_HINT, resolver=which)

    inputs = [str(f) for f in files]
    if _looks_like_fasta(inputs):
        raise RuntimeError(FASTA_INPUT_ERROR)

    base = Path(db_dir) if db_dir else SPECIES_DB_DIR / "D6_sylph"
    dbs = sorted(base.glob("*.syl2db")) if base.is_dir() else []
    if not dbs:
        raise RuntimeError(f"sylph database not found in {base} (download from {_SYLPH_DOCS_URL})")
    db_file = dbs[0]

    with tempfile.TemporaryDirectory(prefix="gside-sylph-") as tmp:
        work = Path(tmp)
        # sylph v1: read sketches ignore -o (genomes only) — they go to a
        # folder given via -d; paired reads must ride the -1/-2 flags.
        sketch_dir = work / "sketches"
        sketch_cmd = [sylph_bin, "sketch", "-t", str(threads)]
        if len(inputs) == 2:
            sketch_cmd += ["-1", inputs[0], "-2", inputs[1]]
        else:
            sketch_cmd += inputs
        sketch_cmd += ["-d", str(sketch_dir)]
        run_checked(sketch_cmd, timeout=3600, name="sylph sketch failed")

        profile = work / "profile.tsv"
        run_checked(
            [sylph_bin, "profile", "-d", str(db_file), str(_find_sketch(work)), "-o", str(profile)],
            timeout=600,
            name="sylph profile failed",
        )
        species = _parse_profile(profile)

    return SylphResult(
        database={"name": "D6_sylph", "path": str(db_file)},
        result={"species": species, "n_species": len(species)},
    )
