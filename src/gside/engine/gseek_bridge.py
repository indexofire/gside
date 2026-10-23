"""Optional gseek bridge for gside's marker BLAST and FASTQ screening.

``GSIDE_ENGINE`` selects the implementation:

- ``legacy`` (default): existing subprocess backends.
- ``gseek``: use gseek's BLAST/minimap2 adapters.
- ``shadow``: run both, log differences, return the legacy result.

The bridge is imported unconditionally (it has no gseek import at module load),
so installations that never set ``GSIDE_ENGINE`` do not need gseek. The
prebuilt-DB paths (panel/skani, mash_refseq, sourmash) are intentionally not
migrated yet: gseek's index model consumes FASTA, while gside ships skani/mash
database directories and a sourmash GTDB zip.
"""

from __future__ import annotations

import importlib
import importlib.util
import logging
import os
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from .hits import Hit

logger = logging.getLogger(__name__)

_VALID_MODES = frozenset({"legacy", "gseek", "shadow"})

__all__ = [
    "aggregate_screen",
    "blast_hits",
    "engine_mode",
    "gseek_available",
    "log_blast_shadow",
    "log_screen_shadow",
    "record_to_hit",
    "screen_reads_gseek",
]


def engine_mode() -> str:
    """Return the configured engine mode (``legacy`` by default)."""
    raw = os.environ.get("GSIDE_ENGINE", "legacy").strip().lower() or "legacy"
    if raw not in _VALID_MODES:
        raise RuntimeError(f"unknown GSIDE_ENGINE={raw!r}; expected one of {sorted(_VALID_MODES)}")
    return raw


def gseek_available() -> bool:
    """True when the optional gseek package can be imported."""
    return importlib.util.find_spec("gseek") is not None


def _require_gseek() -> Any:
    try:
        return importlib.import_module("gseek")
    except ImportError as exc:  # pragma: no cover - depends on optional install
        raise RuntimeError(
            "GSIDE_ENGINE=gseek requires the optional gseek package; "
            "install it with `pip install -e ../gseek` (or `pip install gseek`)"
        ) from exc


def _resolver() -> Any:
    """Use the same pixi/PATH discovery order as gside."""
    from gside.config import pixi_path

    execution = importlib.import_module("gseek.execution")
    ToolResolver = execution.ToolResolver
    paths = [Path(part) for part in pixi_path().split(os.pathsep) if part and Path(part).is_dir()]
    return ToolResolver(extra_paths=paths)


def _coverage(aligned: int, length: int) -> float:
    return 100.0 * aligned / length if length > 0 else 0.0


def record_to_hit(record: Any) -> Hit:
    """Convert a gseek ``AlignmentRecord`` to gside's BLAST ``Hit``.

    gseek coordinates are 0-based half-open; ``Hit`` mirrors BLAST's 1-based
    inclusive coordinates. Coverage deliberately keeps gside's historical
    ``aligned_length / length`` formula rather than gseek's span coverage.
    """
    aligned = int(record.aligned_length or 0)
    query_length = int(record.query_length or 0)
    subject_length = int(record.target_length or 0)
    query_start = int(record.query_start) + 1 if record.query_start is not None else 0
    query_end = int(record.query_end or 0)
    subject_start = int(record.target_start) + 1 if record.target_start is not None else 0
    subject_end = int(record.target_end or 0)

    return Hit(
        query_id=str(record.query_id),
        subject_id=str(record.target_id),
        identity=float(record.identity or 0.0),
        query_coverage=_coverage(aligned, query_length),
        subject_coverage=_coverage(aligned, subject_length),
        evalue=float(record.evalue or 0.0),
        bit_score=float(record.bit_score or 0.0),
        query_start=query_start,
        query_end=query_end,
        subject_start=subject_start,
        subject_end=subject_end,
        strand="-" if bool(getattr(record, "is_reverse", False)) else "+",
        alignment_length=aligned,
        query_length=query_length,
        subject_length=subject_length,
        mismatches=int(record.mismatches or 0),
        backend="blast",
    )


def blast_hits(
    query: Path,
    target_fasta: Path,
    *,
    evalue: float = 1e-5,
    max_targets: int = 500,
    threads: int = 4,
    word_size: int = 11,
    min_identity: float = 0.0,
    min_coverage: float = 0.0,
) -> list[Hit]:
    """Run gseek's BLAST adapter for the marker database.

    ``task="megablast"`` and ``dust=True`` reproduce gside's historical
    ``blastn`` default (BLAST+ defaults blastn to megablast with dust on).
    """
    gseek = _require_gseek()
    result = gseek.align(
        backend="blast",
        query=query,
        target=target_fasta,
        params={
            "task": "megablast",
            "evalue": evalue,
            "max_targets": max_targets,
            "word_size": word_size,
            "dust": True,
        },
        threads=threads,
        resolver=_resolver(),
        index="auto",
    )
    hits = [record_to_hit(record) for record in result.records]
    return [
        hit for hit in hits if hit.identity >= min_identity and hit.query_coverage >= min_coverage
    ]


def _union_length(spans: Iterable[tuple[int, int]]) -> int:
    covered = 0
    reach = -1
    for start, end in sorted(spans):
        if end <= reach:
            continue
        covered += end - max(start, reach)
        reach = end
    return covered


def aggregate_screen(records: Iterable[Any]) -> list[dict[str, Any]]:
    """Aggregate gseek minimap2 records exactly like gside's PAF reducer.

    Breadth is the union of target intervals over gene length; identity is the
    alignment-length-weighted mean of ``nmatch / alen`` (gside PAF semantics),
    not gseek's NM-first identity.
    """
    spans: dict[str, list[tuple[int, int]]] = {}
    lengths: dict[str, int] = {}
    weights: dict[str, list[tuple[float, int]]] = {}

    for record in records:
        aligned = int(record.aligned_length or 0)
        if aligned <= 0:
            continue
        if record.target_start is None or record.target_end is None:
            continue
        target = str(record.target_id)
        spans.setdefault(target, []).append((int(record.target_start), int(record.target_end)))
        lengths[target] = int(record.target_length or 0)
        if record.matches is None:
            continue
        identity = 100.0 * int(record.matches) / aligned
        weights.setdefault(target, []).append((identity, aligned))

    hits: list[dict[str, Any]] = []
    for target, target_spans in spans.items():
        target_length = lengths[target]
        coverage = 100.0 * _union_length(target_spans) / target_length if target_length > 0 else 0.0
        pairs = weights.get(target, [])
        weight_sum = sum(weight for _, weight in pairs)
        identity = (
            sum(value * weight for value, weight in pairs) / weight_sum if weight_sum else 0.0
        )
        parts = target.split("~~~")
        hits.append(
            {
                "gene": parts[1].lower() if len(parts) >= 2 else target.lower(),
                "identity": round(identity, 2),
                "coverage": round(coverage, 1),
                "ref_len": target_length,
                "seqid": target,
            }
        )
    return hits


def screen_reads_gseek(
    reads: str | Path | Sequence[str | Path],
    markers_fasta: str | Path,
    *,
    threads: int = 4,
    preset: str = "sr",
) -> list[dict[str, Any]]:
    """Map FASTQ reads with gseek's minimap2 adapter and aggregate."""
    gseek = _require_gseek()
    files = [Path(reads)] if isinstance(reads, (str, Path)) else [Path(p) for p in reads]
    result = gseek.align(
        backend="minimap2",
        query=tuple(files),
        target=Path(markers_fasta),
        query_kind="fastq",
        paired=len(files) > 1,
        params={"preset": preset, "secondary": "no"},
        threads=threads,
        resolver=_resolver(),
        index="auto",
    )
    return aggregate_screen(result.records)


def log_blast_shadow(legacy: Sequence[Hit], new: Sequence[Hit]) -> None:
    """Log a coarse legacy-vs-gseek BLAST difference for shadow runs."""
    legacy_keys = {(hit.subject_id, round(hit.identity, 2)) for hit in legacy}
    new_keys = {(hit.subject_id, round(hit.identity, 2)) for hit in new}
    if legacy_keys != new_keys:
        logger.warning(
            "gseek shadow BLAST differs: legacy=%d gseek=%d only_legacy=%d only_gseek=%d",
            len(legacy),
            len(new),
            len(legacy_keys - new_keys),
            len(new_keys - legacy_keys),
        )


def log_screen_shadow(legacy: Sequence[dict[str, Any]], new: Sequence[dict[str, Any]]) -> None:
    """Log a coarse legacy-vs-gseek read-screening difference."""
    legacy_genes = {hit["gene"] for hit in legacy}
    new_genes = {hit["gene"] for hit in new}
    if legacy_genes != new_genes:
        logger.warning(
            "gseek shadow read screening differs: legacy_genes=%s gseek_genes=%s",
            sorted(legacy_genes),
            sorted(new_genes),
        )
