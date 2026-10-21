"""Sylph ANI-based metagenomic profiling (gside meta, FASTQ reads only)."""

from __future__ import annotations

import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

from gside.analysis._envelope import METAGENOMIC, MethodEnvelope
from gside.analysis._thresholds import ani_confidence
from gside.config import SPECIES_DB_DIR, which
from gside.engine._env import require_bin, run_checked

_SYLPH_HINT = "sylph not found — pixi add sylph"
FASTA_INPUT_ERROR = "sylph requires FASTQ reads, not contigs — use --mode kraken2 for FASTA input"

_FASTA_SUFFIXES = (".fna", ".fa", ".fasta")


@dataclass
class SylphResult(MethodEnvelope):
    analysis_type: ClassVar[str] = METAGENOMIC
    method: str = "sylph_ani"


def _looks_like_fasta(files: Sequence[str]) -> bool:
    return any(Path(f).suffix in _FASTA_SUFFIXES for f in files)


def _species_from_genome_file(genome_file: str) -> str:
    """`Escherichia_coli.fna` → `Escherichia coli`; GCF accessions resolved via metadata."""
    stem = Path(genome_file).stem
    if stem.startswith("GCF_") or stem.startswith("GCA_"):
        return stem
    return stem.replace("_", " ")


def _load_taxa_map(db_dir: Path) -> dict[str, str]:
    """Load filename→species from metadata.tsv if present near the sylph DB."""
    for candidate in (
        db_dir / "metadata.tsv",
        db_dir.parent / "metadata.tsv",
        db_dir.parent / "D2_ani" / "metadata.tsv",
    ):
        if candidate.is_file():
            mapping: dict[str, str] = {}
            for ln in candidate.read_text().splitlines()[1:]:
                cols = ln.split("\t")
                if len(cols) >= 2:
                    full = cols[0].strip()
                    species = cols[1].strip()
                    mapping[full] = species
                    mapping[Path(full).stem] = species
            return mapping
    return {}


def _find_sketch(work: Path) -> Path:
    """Locate the read sketch sylph wrote (v0 emits .syldata, v1 emits .sylsp
    inside the ``-d`` sketches directory)."""
    candidates = sorted(work.rglob("*.sylsp")) + sorted(work.rglob("*.syldata"))
    if not candidates:
        raise RuntimeError("sylph sketch produced no .syldata/.sylsp output")
    return candidates[0]


def _parse_profile(path: Path, taxa_map: dict[str, str] | None = None) -> list[dict[str, Any]]:
    """sylph profile TSV → species rows sorted by taxonomic abundance descending."""
    mapping = taxa_map or {}
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    species: list[dict[str, Any]] = []
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t"), strict=False))
        try:
            ani = float(row.get("Adjusted_ANI", 0.0))
            abundance = float(row.get("Taxonomic_abundance", 0.0)) / 100.0
            coverage = float(row.get("Eff_cov", 0.0))
        except ValueError:
            continue
        raw_name = row.get("Genome_file", "")
        stem = Path(raw_name).stem
        name = mapping.get(raw_name, mapping.get(stem, _species_from_genome_file(raw_name)))
        species.append(
            {
                "name": name,
                "ani": ani,
                "abundance": abundance,
                "coverage": coverage,
                "confidence": ani_confidence(ani),
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
    dbs: list[Path] = []
    if base.is_dir():
        # .syl2db is the v1 two-stage format and is preferred by sylph itself.
        dbs = sorted(base.glob("*.syl2db")) + sorted(base.glob("*.syldb"))
    if not dbs:
        raise RuntimeError(
            f"sylph database not found in {base} (run: gside db setup --tier meta_sylph)"
        )
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
        taxa_map = _load_taxa_map(base)
        species = _parse_profile(profile, taxa_map)

    return SylphResult(
        database={"name": "D6_sylph", "path": str(db_file)},
        result={"species": species, "n_species": len(species)},
    )
