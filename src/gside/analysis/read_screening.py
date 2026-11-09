"""FASTQ marker screening via minimap2 (gspec reads-mode recipe).

Reads are mapped against the marker FASTA and per-gene breadth is the union
of reference intervals across all primary alignments divided by gene length
(the same aggregate-breadth math as gspec's reads mode). Identity is the
alignment-length-weighted mean of ``nmatch / alen`` from the PAF rows. The
multigene rule engine then runs unchanged, so a species needs the same
per-gene identity/breadth as it would from an assembly.

Preset: short reads (mean < 300 bp, averaged over mates) use ``-x sr``; longer reads use
``-x map-ont``. ``--secondary=no`` keeps only the primary alignment of each
read; identity/breadth thresholds live in the rule engine, not here.
"""

from __future__ import annotations

import gzip
from collections.abc import Sequence
from pathlib import Path
from typing import Any, TextIO

from gside.config import which
from gside.engine._env import require_bin, run_checked
from gside.engine.gseek_bridge import engine, log_screen_shadow, screen_reads_gseek

_FASTQ_SUFFIXES = (".fq", ".fastq")
_SHORT_READ_CUTOFF = 300  # bp: below -> sr preset, at/above -> map-ont
_SAMPLE_READS = 500


def is_fastq(path: str | Path) -> bool:
    """True for ``.fq``/``.fastq`` (optionally ``.gz``)."""
    name = str(path).lower()
    for suffix in (".gz", ".bz2"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name.endswith(_FASTQ_SUFFIXES)


def _open_text(path: str | Path) -> TextIO:
    return gzip.open(path, "rt") if str(path).lower().endswith(".gz") else open(path)


def mean_read_length(path: str | Path, sample: int = _SAMPLE_READS) -> float:
    """Mean length of the first ``sample`` reads (0.0 for an empty FASTQ)."""
    total = 0
    count = 0
    with _open_text(path) as fh:
        while count < sample:
            header = fh.readline()
            if not header:
                break
            if not header.startswith("@"):
                continue
            seq = fh.readline().strip()
            fh.readline()  # plus line
            fh.readline()  # quality
            if seq:
                total += len(seq)
                count += 1
    return total / count if count else 0.0


def _parse_paf_line(line: str) -> dict[str, Any] | None:
    """One primary PAF row → dict; None for malformed/secondary rows."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    tags: dict[str, str] = {}
    for raw in fields[12:]:
        key, _, value = raw.partition(":")
        tags[key] = value
    if tags.get("tp", "A:P").endswith(":S"):
        return None
    try:
        return {
            "qname": fields[0],
            "tname": fields[5],
            "tlen": int(fields[6]),
            "tstart": int(fields[7]),
            "tend": int(fields[8]),
            "nmatch": int(fields[9]),
            "alen": int(fields[10]),
            "mapq": int(fields[11]),
        }
    except ValueError:
        return None


def _union_length(spans: list[tuple[int, int]]) -> int:
    """Total length of the union of half-open [start, end) intervals."""
    covered = 0
    reach = -1
    for start, end in sorted(spans):
        if end <= reach:
            continue
        covered += end - max(start, reach)
        reach = end
    return covered


def _aggregate(rows: list[dict[str, Any]], min_mapq: int = 0) -> list[dict[str, Any]]:
    """Per-target breadth (union), alignment-weighted identity, and mean depth.

    ``mean_depth`` follows gspec's reads-mode contract: the summed length of
    all primary target intervals divided by the target length.  It is not a
    read count, so genes of different lengths are directly comparable.
    """
    spans: dict[str, list[tuple[int, int]]] = {}
    lengths: dict[str, int] = {}
    weights: dict[str, list[tuple[float, int]]] = {}
    qnames: dict[str, set[str]] = {}
    mapqs: dict[str, list[int]] = {}

    for row in rows:
        if row["alen"] <= 0 or row.get("mapq", 0) < min_mapq:
            continue
        target = row["tname"]
        spans.setdefault(target, []).append((row["tstart"], row["tend"]))
        lengths[target] = row["tlen"]
        qnames.setdefault(target, set()).add(row["qname"])
        mapqs.setdefault(target, []).append(row.get("mapq", 0))
        identity = 100.0 * row["nmatch"] / row["alen"]
        weights.setdefault(target, []).append((identity, row["alen"]))

    hits: list[dict[str, Any]] = []
    for target, target_spans in spans.items():
        tlen = lengths[target]
        covered = _union_length(target_spans)
        breadth = 100.0 * covered / tlen if tlen else 0.0
        aligned_bases = sum(end - start for start, end in target_spans)
        mean_depth = aligned_bases / tlen if tlen else 0.0
        pairs = weights[target]
        weight_sum = sum(weight for _, weight in pairs)
        identity = sum(value * weight for value, weight in pairs) / weight_sum
        target_mapqs = mapqs[target]
        from gside.marker_db import gene_from_seqid

        parts = [target, gene_from_seqid(target), ""]
        hits.append(
            {
                "gene": parts[1].lower() if len(parts) >= 2 else target.lower(),
                "identity": round(identity, 2),
                "coverage": round(breadth, 1),
                "mean_depth": round(mean_depth, 2),
                "reads": len(qnames[target]),
                "aligned_bases": aligned_bases,
                "mapq_mean": round(sum(target_mapqs) / len(target_mapqs), 1),
                "ref_len": tlen,
                "seqid": target,
            }
        )
    return hits


def _screen_reads_legacy(
    files: Sequence[str | Path],
    markers_fasta: str | Path,
    threads: int,
    preset: str,
    min_mapq: int = 20,
) -> list[dict[str, Any]]:
    minimap2 = require_bin(
        "minimap2",
        hint="minimap2 not found — install via pixi (conda) or point GSIDE_PIXI_BIN at a bin dir",
        resolver=which,
    )
    cmd = [
        minimap2,
        "-x",
        preset,
        "-c",
        "-t",
        str(threads),
        "--secondary=no",
        str(markers_fasta),
        *[str(f) for f in files],
    ]
    proc = run_checked(cmd, timeout=3600, name="minimap2 marker screening failed")
    rows = [row for line in proc.stdout.splitlines() if (row := _parse_paf_line(line))]
    return _aggregate(rows, min_mapq=min_mapq)


def screen_reads(
    reads: str | Path | Sequence[str | Path],
    markers_fasta: str | Path,
    threads: int = 4,
    min_mapq: int = 20,
) -> list[dict[str, Any]]:
    """Map one FASTQ file (or both mates) to the marker FASTA.

    Paired mates are processed in one minimap2 invocation, so per-gene
    breadth is the union of intervals from both files.  ``-c`` is required:
    without base-level alignment minimap2's PAF columns 10/11 are approximate,
    and they must not be treated as exact matches/alignment length.

    ``GSIDE_ENGINE`` switches between gseek's minimap2 adapter (default),
    the legacy subprocess path, and a shadow comparison.
    """
    files = [reads] if isinstance(reads, (str, Path)) else list(reads)
    lengths = [mean_read_length(f) for f in files] or [0.0]
    mean = sum(lengths) / len(lengths)
    preset = "sr" if mean < _SHORT_READ_CUTOFF else "map-ont"

    def legacy() -> list[dict[str, Any]]:
        return _screen_reads_legacy(files, markers_fasta, threads, preset, min_mapq=min_mapq)

    def compare(old: list[dict[str, Any]], new: list[dict[str, Any]]) -> bool:
        log_screen_shadow(old, new)
        return True

    return engine().run(
        legacy,
        lambda: screen_reads_gseek(
            files, markers_fasta, threads=threads, preset=preset, min_mapq=min_mapq
        ),
        compare=compare,
    )


__all__ = ["is_fastq", "mean_read_length", "screen_reads"]
