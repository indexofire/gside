"""Optional gseek bridge for gside's marker BLAST and FASTQ screening.

``GSIDE_ENGINE`` selects the implementation:

- ``legacy`` (default): existing subprocess backends.
- ``gseek``: use gseek's BLAST/minimap2 adapters.
- ``shadow``: run both, log differences, return the legacy result.

The bridge is imported unconditionally (it has no gseek import at module load),
so installations that never set ``GSIDE_ENGINE`` do not need gseek. Marker
BLAST/FASTQ screening, mash distances (assembly and reads), sourmash gather, and
skani panel searches use gseek. gseek external references accept pre-built
files and database directories (sourmash rocksdb, skani sketch DB).
"""

from __future__ import annotations

import functools
import importlib
import importlib.util
import logging
import os
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeVar

from .hits import Hit

if TYPE_CHECKING:
    from gseek.compat import Engine

logger = logging.getLogger(__name__)

_VALID_MODES = frozenset({"legacy", "gseek", "shadow"})

__all__ = [
    "aggregate_screen",
    "blast_hits",
    "engine",
    "engine_mode",
    "gseek_available",
    "log_blast_shadow",
    "log_mash_shadow",
    "log_skani_shadow",
    "log_screen_shadow",
    "log_sourmash_shadow",
    "log_sylph_shadow",
    "mash_distance",
    "mash_distance_reads",
    "record_to_hit",
    "screen_reads_gseek",
    "skani_search",
    "sourmash_gather",
    "sylph_profile",
]


T = TypeVar("T")


class _LegacyMode:
    value = "legacy"


class _LegacyEngine:
    """No-gseek fallback; legacy installs never import the optional package."""

    mode = _LegacyMode()

    def run(
        self,
        legacy: Callable[[], T],
        gseek: Callable[[], T],
        *,
        compare: Callable[[T, T], bool] | None = None,
    ) -> T:
        return legacy()


def engine() -> _LegacyEngine | Engine:
    """Return the shared switch, importing gseek only for gseek/shadow modes."""
    raw = os.environ.get("GSIDE_ENGINE", "legacy").strip().lower() or "legacy"
    if raw not in _VALID_MODES:
        raise RuntimeError(f"unknown GSIDE_ENGINE={raw!r}; expected one of {sorted(_VALID_MODES)}")
    if raw == "legacy":
        return _LegacyEngine()
    try:
        from gseek.compat import Engine
    except ImportError as exc:
        raise RuntimeError("GSIDE_ENGINE=gseek requires the optional gseek package") from exc
    return Engine.from_env("GSIDE_ENGINE", logger=logger)


def engine_mode() -> str:
    """Return the configured engine mode (``legacy`` by default)."""
    return engine().mode.value


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


@functools.lru_cache(maxsize=1)
def _resolver() -> Any:
    """Use the same pixi/PATH discovery order as gside.

    Tool version probing is cached by the resolver instance; reusing it avoids
    one ``--version`` subprocess per alignment call in batch runs.
    """
    from gside.config import pixi_path

    execution = importlib.import_module("gseek.execution")
    ToolResolver = execution.ToolResolver
    paths = [Path(part) for part in pixi_path().split(os.pathsep) if part and Path(part).is_dir()]
    return ToolResolver(extra_paths=paths)


_REFERENCE_DIGESTS: dict[str, tuple[int, str]] = {}


def _reference_digest(path: Path) -> tuple[int, str]:
    """Return (size, sha256) for an external reference, cached per process."""
    provenance = importlib.import_module("gseek.provenance")
    directory_size = provenance.directory_size
    sha256_file = provenance.sha256_file
    sha256_tree = provenance.sha256_tree

    resolved = path.resolve()
    key = str(resolved)
    cached = _REFERENCE_DIGESTS.get(key)
    if cached is not None:
        return cached
    if resolved.is_dir():
        value = (directory_size(resolved), sha256_tree(resolved))
    else:
        value = (resolved.stat().st_size, sha256_file(resolved))
    _REFERENCE_DIGESTS[key] = value
    return value


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
        record_manifest=False,
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
        record_manifest=False,
    )
    return aggregate_screen(result.records)


def _mash_records(gseek: Any, **kwargs: Any) -> list[Any]:
    from gside.engine.backends.kmer import KmerDistance

    result = gseek.distance(
        backend="mash",
        external_reference=True,
        external_sha256=kwargs.pop("external_sha256", None),
        external_size=kwargs.pop("external_size", None),
        resolver=_resolver(),
        index="none",
        record_manifest=False,
        **kwargs,
    )
    return [
        KmerDistance(
            query_id=record.query_id,
            reference_id=record.reference_id,
            distance=float(record.distance or 0.0),
            pvalue=float(record.pvalue or 0.0),
            shared_hashes=int(record.shared_hashes or 0),
            total_hashes=int(record.total_hashes or 0),
            backend="mash",
        )
        for record in result.records
    ]


def mash_distance(
    query: str | Path,
    reference_msh: str | Path,
    *,
    threads: int = 4,
    max_distance: float = 1.0,
) -> list[Any]:
    """Run gseek mash against a pre-built ``.msh`` reference.

    The query is FASTA/FASTQ; gseek sketches it with the default Mash settings
    and passes the external reference straight to ``mash dist``.
    """
    reference = Path(reference_msh)
    size, digest = _reference_digest(reference)
    return _mash_records(
        _require_gseek(),
        query=Path(query),
        reference=reference,
        external_sha256=digest,
        external_size=size,
        params={"max_distance": max_distance},
        threads=threads,
    )


def mash_distance_reads(
    reads: str | Path | Sequence[str | Path],
    reference_msh: str | Path,
    *,
    threads: int = 4,
    max_distance: float = 1.0,
) -> list[Any]:
    """Run gseek mash on FASTQ reads with gside's ``-m 2 -r`` sketch flags."""
    files = [Path(reads)] if isinstance(reads, (str, Path)) else [Path(p) for p in reads]
    reference = Path(reference_msh)
    size, digest = _reference_digest(reference)
    return _mash_records(
        _require_gseek(),
        query=tuple(files),
        reference=reference,
        external_sha256=digest,
        external_size=size,
        params={"max_distance": max_distance, "reads": True},
        threads=threads,
    )


def sourmash_gather(
    contigs: str | Path | Sequence[str | Path],
    reference_db: str | Path,
    *,
    threads: int = 4,
) -> list[dict[str, Any]]:
    """Run gseek sourmash gather against a pre-built GTDB zip/signature.

    Returns rows shaped like the legacy gather CSV (``name``,
    ``f_unique_weighted``, ``average_containment_ani``) so gside's lineage
    mapping stays unchanged.
    """
    gseek = _require_gseek()
    files = [Path(contigs)] if isinstance(contigs, (str, Path)) else [Path(p) for p in contigs]
    reference = Path(reference_db)
    size, digest = _reference_digest(reference)
    rows: list[dict[str, Any]] = []
    # sourmash gather requires exactly one query signature, so batch multi-file
    # samples as one gseek run per file and merge the gather rows.
    for path in files:
        try:
            result = gseek.distance(
                backend="sourmash",
                query=path,
                reference=reference,
                external_reference=True,
                external_sha256=digest,
                external_size=size,
                params={
                    "mode": "gather",
                    "k": 31,
                    "scaled": 1000,
                    # gside relies on sourmash's CLI default; gseek's default
                    # is 50 bp, which admits extra low-weight partitions.
                    "threshold_bp": 50000,
                },
                threads=threads,
                resolver=_resolver(),
                index="none",
                record_manifest=False,
            )
        except Exception as exc:  # gseek errors are typed; gside's caller expects RuntimeError
            raise RuntimeError(f"gseek sourmash gather failed for {path}: {exc}") from exc
        for record in result.records:
            name = record.tags.get("reference_name") or record.reference_id
            ani = (record.ani or 0.0) / 100.0
            rows.append(
                {
                    "name": name,
                    "f_unique_weighted": record.tags.get("f_unique_weighted", "0"),
                    "average_containment_ani": str(ani),
                }
            )
    return rows


def _mash_reference_id(item: Any) -> str:
    if hasattr(item, "reference_id"):
        return str(item.reference_id)
    return str(item[0])


def _mash_distance(item: Any) -> float:
    if hasattr(item, "distance"):
        return float(item.distance)
    # legacy tuple is (reference_id, identity); identity = 1 - distance
    return 1.0 - float(item[1])


def skani_search(
    query: str | Path,
    reference_db: str | Path,
    *,
    threads: int = 4,
) -> list[Any]:
    """Run gseek skani search against a pre-built database directory.

    gside historically used the query-side aligned fraction
    (``Align_fraction_query``) and the ``Ref_file`` basename. gseek's normalized
    record carries the reference fraction in ``aligned_fraction``, so map the
    raw query fraction tag back to gside's ``AniHit``.
    """
    from gside.engine.backends.skani import AniHit

    gseek = _require_gseek()
    reference = Path(reference_db)
    size, digest = _reference_digest(reference)
    result = gseek.distance(
        backend="skani",
        query=Path(query),
        reference=reference,
        external_reference=True,
        external_sha256=digest,
        external_size=size,
        threads=threads,
        resolver=_resolver(),
        index="none",
        record_manifest=False,
    )
    hits: list[AniHit] = []
    for record in result.records:
        ref = record.tags.get("reference_file") or record.reference_id
        aligned = float(record.tags.get("query_fraction", record.aligned_fraction or 0.0))
        hits.append(
            AniHit(
                ref=ref,
                ani=float(record.ani or 0.0),
                aligned_fraction=aligned,
            )
        )
    return hits


def log_sylph_shadow(legacy: Sequence[Any], new: Sequence[Any]) -> None:
    def key(rows: Sequence[Any]) -> list[tuple[str, float, float, float]]:
        return sorted(
            (
                str(row.get("name", row.get("genome_file", ""))),
                round(float(row.get("ani", 0.0)), 3),
                round(float(row.get("abundance", 0.0)), 4),
                round(float(row.get("coverage", 0.0)), 3),
            )
            for row in rows
        )

    legacy_keys, new_keys = key(legacy), key(new)
    if legacy_keys != new_keys:
        logger.warning(
            "gseek shadow sylph differs: legacy=%d gseek=%d only_legacy=%d only_gseek=%d",
            len(legacy_keys),
            len(new_keys),
            len(set(legacy_keys) - set(new_keys)),
            len(set(new_keys) - set(legacy_keys)),
        )


def sylph_profile(
    reads: str | Path | Sequence[str | Path],
    reference_db: str | Path,
    *,
    threads: int = 4,
) -> list[dict[str, Any]]:
    """Run gseek's sylph backend against a pre-built ``.syldb``/``.syl2db``.

    Returns rows shaped like the legacy profile parser *before* species-name
    mapping (``genome_file``, ``ani``, ``abundance``, ``coverage``): the
    metadata mapping and ANI-confidence thresholds stay in gside.
    """
    gseek = _require_gseek()
    files = [Path(reads)] if isinstance(reads, (str, Path)) else [Path(p) for p in reads]
    reference = Path(reference_db)
    size, digest = _reference_digest(reference)
    result = gseek.distance(
        backend="sylph",
        query=files,
        reference=reference,
        external_reference=True,
        external_sha256=digest,
        external_size=size,
        threads=threads,
        resolver=_resolver(),
        index="none",
        record_manifest=False,
    )
    rows: list[dict[str, Any]] = []
    for record in result.records:
        tags = record.tags
        abundance_raw = tags.get("taxonomic_abundance", "0") or "0"
        try:
            abundance = float(abundance_raw) / 100.0
        except ValueError:
            abundance = 0.0
        rows.append(
            {
                "genome_file": str(tags.get("genome_file", "")),
                "ani": float(record.ani or 0.0),
                "abundance": abundance,
                # gseek stores Eff_cov as a 0..1 fraction; gside keeps percent.
                "coverage": float(record.aligned_fraction or 0.0) * 100.0,
            }
        )
    rows.sort(key=lambda row: row["abundance"], reverse=True)
    return rows


def log_skani_shadow(legacy: Sequence[Any], new: Sequence[Any]) -> None:
    legacy_keys = sorted(
        (hit.ref, round(hit.ani, 3), round(hit.aligned_fraction, 3)) for hit in legacy
    )
    new_keys = sorted((hit.ref, round(hit.ani, 3), round(hit.aligned_fraction, 3)) for hit in new)
    if legacy_keys != new_keys:
        logger.warning(
            "gseek shadow skani differs: legacy=%d gseek=%d only_legacy=%d only_gseek=%d",
            len(legacy_keys),
            len(new_keys),
            len(set(legacy_keys) - set(new_keys)),
            len(set(new_keys) - set(legacy_keys)),
        )


def log_mash_shadow(legacy: Sequence[Any], new: Sequence[Any]) -> None:
    legacy_top = sorted(legacy, key=_mash_distance)[:10]
    new_top = sorted(new, key=_mash_distance)[:10]
    if [_mash_reference_id(item) for item in legacy_top] != [
        _mash_reference_id(item) for item in new_top
    ]:
        logger.warning(
            "gseek shadow mash differs: legacy_top=%s gseek_top=%s",
            [_mash_reference_id(item) for item in legacy_top],
            [_mash_reference_id(item) for item in new_top],
        )


def log_sourmash_shadow(legacy: Sequence[dict[str, Any]], new: Sequence[dict[str, Any]]) -> None:
    legacy_names = sorted(row.get("name", "") for row in legacy)
    new_names = sorted(row.get("name", "") for row in new)
    if legacy_names != new_names:
        logger.warning(
            "gseek shadow sourmash differs: legacy_top=%s gseek_top=%s",
            legacy_names,
            new_names,
        )


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
