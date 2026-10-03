"""ANI-based species identification (species-id plan A: panel / skani_gtdb / mash_refseq)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from gside.config import SPECIES_DB_DIR
from gside.engine.backends.skani import AniHit

_DB_BY_MODE = {
    "panel": "D2_ani",
    "skani_gtdb": "D2_skani_gtdb",
    "mash_refseq": "D3_mash",
}

_ANI_HIGH = 95.0
_ANI_MEDIUM = 93.0
_AF_MIN = 0.65
_MASH_IDENTITY_HIGH = 0.97
_TOP_N = 10


@dataclass
class AniIdResult:
    method: str
    database: dict[str, str]
    result: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "analysis_type": "species_identification",
            "method": self.method,
            "database": self.database,
            "result": self.result,
        }


def _hits_from_rows(rows: list[tuple[str, float, float]]) -> list[AniHit]:
    return [AniHit(ref=r, ani=a, aligned_fraction=af) for r, a, af in rows]


def _skani_search(query: str | Path, db: Path) -> list[AniHit]:
    from gside.engine.backends.skani import SkaniBackend

    return SkaniBackend().search(Path(query), db)


def _mash_dist(query: str | Path, msh: Path) -> list[tuple[str, float]]:
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
    results = MashBackend(binary=mash).distance(Path(msh), Path(query), max_distance=1.0)
    return [(r.reference_id, 1.0 - r.distance) for r in results]


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


def identify_by_ani(
    contigs: str | Path,
    mode: str,
    db_dir: str | Path | None = None,
) -> AniIdResult:
    if mode not in _DB_BY_MODE:
        raise ValueError(f"unsupported ANI mode {mode!r}; expected one of {sorted(_DB_BY_MODE)}")
    db_name = _DB_BY_MODE[mode]
    base = Path(db_dir) if db_dir else SPECIES_DB_DIR / db_name
    taxa = _load_taxa_map(base)
    database = {"name": db_name, "version": _db_version(base, db_name)}

    if mode == "mash_refseq":
        msh = base / "mash.msh"
        if not msh.exists():
            msh = base / "payload.bin"
        if not msh.exists():
            raise RuntimeError(
                f"mash sketch DB not found in {base} (run: gside db setup --tier mash)"
            )
        mash_rows = _mash_dist(contigs, msh)
        mash_top = sorted(mash_rows, key=lambda r: -r[1])[:_TOP_N]
        best_ref, best_identity = mash_top[0] if mash_top else ("", 0.0)
        species = taxa.get(best_ref, _mash_species(best_ref)) if best_ref else "Unknown"
        confidence = "high" if best_identity >= _MASH_IDENTITY_HIGH else "low"
        result = {
            "species": species or "Unknown",
            "confidence": confidence,
            "identity": best_identity,
            "top_hits": [{"genome": r, "identity": i} for r, i in mash_top],
        }
        if 0.90 <= best_identity < _MASH_IDENTITY_HIGH:
            result["confidence"] = "medium"
        return AniIdResult(method=mode, database=database, result=result)

    skani_db = base / "panel.sketch" if (base / "panel.sketch").exists() else base
    skani_hits = sorted(
        _skani_search(contigs, skani_db), key=lambda h: (-h.ani, -h.aligned_fraction)
    )
    top = skani_hits[:_TOP_N]
    best = top[0] if top else None

    species, confidence = "Unknown", "low"
    if best is not None:
        if best.ani >= _ANI_HIGH and best.aligned_fraction >= _AF_MIN:
            species, confidence = taxa.get(best.ref, best.ref), "high"
        elif _ANI_MEDIUM <= best.ani < _ANI_HIGH and best.aligned_fraction >= _AF_MIN:
            species, confidence = taxa.get(best.ref, best.ref), "medium"

    result = {
        "species": species,
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
