"""ANI-based species identification (species-id plan A: panel / mash_refseq)."""

from __future__ import annotations

import json
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from gside.analysis._envelope import MethodEnvelope
from gside.analysis._thresholds import ANI_HIGH, ANI_MEDIUM
from gside.analysis.read_screening import is_fastq
from gside.analysis.species_name import apply_species_name
from gside.config import SPECIES_DB_DIR
from gside.engine.backends.skani import AniHit

_DB_BY_MODE = {
    "panel": "D2_ani",
    "mash_refseq": "D3_mash",
}

_AF_MIN = 0.65
_MASH_IDENTITY_HIGH = 0.97
_MASH_IDENTITY_MEDIUM = 0.90
_TOP_N = 10


@dataclass
class AniIdResult(MethodEnvelope):
    pass


def _hits_from_rows(rows: list[tuple[str, float, float]]) -> list[AniHit]:
    return [AniHit(ref=r, ani=a, aligned_fraction=af) for r, a, af in rows]


def _skani_search(query: str | Path, db: Path, threads: int = 4) -> list[AniHit]:
    from gside.engine.backends.skani import SkaniBackend

    return SkaniBackend(threads=threads).search(Path(query), db)


def _mash_dist_legacy(query: str | Path, msh: Path, threads: int) -> list[tuple[str, float]]:
    from gside.config import which
    from gside.engine._env import require_bin
    from gside.engine.backends.kmer import MashBackend

    mash = require_bin(
        "mash",
        hint="mash not found — install via pixi (conda) or point GSIDE_PIXI_BIN at a bin dir",
        resolver=which,
    )
    # mash CLI is `dist <reference> <query>`: output column 0 carries the
    # reference (sketch DB genome) IDs, so the sketch goes in the first slot.
    # MashBackend.distance() places its first positional there; max_distance=1.0
    # is mash's own default cutoff, i.e. no extra filtering vs. a bare `mash dist`.
    results = MashBackend(binary=mash, threads=threads).distance(
        Path(msh), Path(query), max_distance=1.0
    )
    return [(r.reference_id, 1.0 - r.distance) for r in results]


def _mash_dist(query: str | Path, msh: Path, threads: int = 4) -> list[tuple[str, float]]:
    from gside.engine.gseek_bridge import (
        engine_mode,
        log_mash_shadow,
        mash_distance,
    )

    mode = engine_mode()
    # Read sketches are still produced by the legacy ``-m 2 -r`` path; only a
    # sequence query (assembly FASTA) is handed to gseek.
    if mode in {"gseek", "shadow"} and Path(query).suffix.lower() != ".msh":
        new = mash_distance(query, msh, threads=threads, max_distance=1.0)
        if mode == "gseek":
            return [(hit.reference_id, 1.0 - hit.distance) for hit in new]
        legacy = _mash_dist_legacy(query, msh, threads)
        log_mash_shadow(legacy, new)
        return legacy
    return _mash_dist_legacy(query, msh, threads)


def _mash_dist_reads_legacy(
    reads: str | Path | Sequence[str | Path], msh: Path, threads: int
) -> list[tuple[str, float]]:
    """Sketch FASTQ reads (error-filtered) and compare against the mash DB.

    ``-m 2`` drops single-copy k-mers (sequencing errors) and ``-r`` estimates
    genome size from k-mer content so p-values stay meaningful. Mash then
    compares a read sketch to the reference sketch: valid for species-level
    calls at >=10x coverage, but strain-level hits are unreliable and the
    distance inflates with coverage/error rate (Mash paper; mash docs).
    """
    from gside.config import which
    from gside.engine._env import require_bin, run_checked

    mash = require_bin(
        "mash",
        hint="mash not found — install via pixi (conda) or point GSIDE_PIXI_BIN at a bin dir",
        resolver=which,
    )
    files = [reads] if isinstance(reads, (str, Path)) else list(reads)
    with tempfile.TemporaryDirectory(prefix="gside-mash-") as tmp:
        sketch = Path(tmp) / "reads.msh"
        run_checked(
            [mash, "sketch", "-m", "2", "-r", "-o", str(sketch), *[str(f) for f in files]],
            timeout=3600,
            name="mash sketch failed",
        )
        return _mash_dist(sketch, msh, threads)


def _mash_dist_reads(
    reads: str | Path | Sequence[str | Path], msh: Path, threads: int = 4
) -> list[tuple[str, float]]:
    from gside.engine.gseek_bridge import (
        engine_mode,
        log_mash_shadow,
        mash_distance_reads,
    )

    mode = engine_mode()
    if mode in {"gseek", "shadow"}:
        new = mash_distance_reads(reads, msh, threads=threads, max_distance=1.0)
        if mode == "gseek":
            return [(hit.reference_id, 1.0 - hit.distance) for hit in new]
        legacy = _mash_dist_reads_legacy(reads, msh, threads)
        log_mash_shadow(legacy, new)
        return legacy
    return _mash_dist_reads_legacy(reads, msh, threads)


def _mash_species(ref_id: str) -> str:
    return ref_id


def _load_taxa_map(db_dir: Path) -> dict[str, str]:
    for candidate in ("metadata.tsv", "taxa_map.tsv"):
        path = db_dir / candidate
        if path.is_file():
            mapping: dict[str, str] = {}
            for ln in path.read_text(encoding="utf-8").splitlines()[1:]:
                cols = ln.split("\t")
                if len(cols) >= 2:
                    mapping[cols[0].strip()] = cols[1].strip()
            return mapping
    return {}


def _db_version(db_dir: Path, name: str) -> str:
    manifest = db_dir.parent / "manifests" / f"{name}.json"
    if manifest.is_file():
        try:
            checksum = json.loads(manifest.read_text()).get("checksum", "")
            return checksum[:8] if checksum else "unknown"
        except json.JSONDecodeError:
            return "unknown"
    return "unknown"


def _mash_result(
    base: Path,
    contigs: str | Path | Sequence[str | Path],
    taxa: dict[str, str],
    threads: int,
) -> dict[str, Any]:
    """mash_refseq: nearest RefSeq sketch by MinHash identity."""
    msh = base / "mash.msh"
    if not msh.exists():
        msh = base / "payload.bin"
    if not msh.exists():
        raise RuntimeError(f"mash sketch DB not found in {base} (run: gside db setup --tier mash)")

    files = [contigs] if isinstance(contigs, (str, Path)) else list(contigs)
    if is_fastq(files[0]):
        mash_rows = _mash_dist_reads(files, msh, threads)
    else:
        mash_rows = _mash_dist(files[0], msh, threads)
    mash_top = sorted(mash_rows, key=lambda r: -r[1])[:_TOP_N]
    best_ref, best_identity = mash_top[0] if mash_top else ("", 0.0)
    raw_species = taxa.get(best_ref, _mash_species(best_ref)) if best_ref else "Unknown"

    confidence = "high" if best_identity >= _MASH_IDENTITY_HIGH else "low"
    if _MASH_IDENTITY_MEDIUM <= best_identity < _MASH_IDENTITY_HIGH:
        confidence = "medium"
    result = {
        "confidence": confidence,
        "identity": best_identity,
        "top_hits": [{"genome": r, "identity": i} for r, i in mash_top],
    }
    return apply_species_name(result, raw_species or "Unknown")


def _panel_result(
    base: Path,
    contigs: str | Path | Sequence[str | Path],
    taxa: dict[str, str],
    threads: int,
) -> dict[str, Any]:
    """panel: skani ANI against the curated reference panel."""
    query = contigs if isinstance(contigs, (str, Path)) else contigs[0]
    skani_db = base / "panel.sketch" if (base / "panel.sketch").exists() else base
    skani_hits = sorted(
        _skani_search(query, skani_db, threads), key=lambda h: (-h.ani, -h.aligned_fraction)
    )
    top = skani_hits[:_TOP_N]
    best = top[0] if top else None

    species, confidence = "Unknown", "low"
    if best is not None and best.aligned_fraction >= _AF_MIN:
        if best.ani >= ANI_HIGH:
            species, confidence = taxa.get(best.ref, best.ref), "high"
        elif best.ani >= ANI_MEDIUM:
            species, confidence = taxa.get(best.ref, best.ref), "medium"

    result = {
        "confidence": confidence,
        "ani": best.ani if best else None,
        "aligned_fraction": best.aligned_fraction if best else None,
        "top_hits": [
            {
                "genome": h.ref,
                "species": taxa.get(h.ref, h.ref),
                "ani": h.ani,
                "af": h.aligned_fraction,
            }
            for h in top
        ],
    }
    return apply_species_name(result, species)


def identify_by_ani(
    contigs: str | Path | Sequence[str | Path],
    mode: str,
    db_dir: str | Path | None = None,
    threads: int = 4,
) -> AniIdResult:
    if mode not in _DB_BY_MODE:
        raise ValueError(f"unsupported ANI mode {mode!r}; expected one of {sorted(_DB_BY_MODE)}")
    db_name = _DB_BY_MODE[mode]
    base = Path(db_dir) if db_dir else SPECIES_DB_DIR / db_name
    taxa = _load_taxa_map(base)
    database = {"name": db_name, "version": _db_version(base, db_name)}

    if mode == "mash_refseq":
        result = _mash_result(base, contigs, taxa, threads)
    else:
        result = _panel_result(base, contigs, taxa, threads)
    return AniIdResult(method=mode, database=database, result=result)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="ANI species identifier")
    parser.add_argument("contigs")
    parser.add_argument("--mode", required=True, choices=sorted(_DB_BY_MODE))
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    res = identify_by_ani(args.contigs, args.mode)
    if args.json:
        print(json.dumps(res.to_dict(), ensure_ascii=False, indent=2))
    else:
        r = res.result
        print(f"Species: {r['species']} ({r['confidence']}) via {res.method}")
        for h in r.get("top_hits", [])[:5]:
            print(f"  {h['genome']}: {h.get('ani', h.get('identity'))}")


if __name__ == "__main__":
    main()
