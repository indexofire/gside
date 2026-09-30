"""gside — Genome Species IDentification Engine.

Multi-layer species identification from contigs:
  marker       multigene combination rules (blastn against markers_v2)
  panel        ANI against the curated reference panel (skani)
  mash_refseq  MinHash distance against RefSeq sketch (mash)
  sourmash     sourmash GTDB gather
  all          run available methods and arbitrate (ANI layer > marker layer)

Output: single JSON verdict per method (GOM-compatible contract).
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any


def _run_species(args: argparse.Namespace) -> int:
    contigs = args.contigs
    results: dict[str, Any] = {}

    modes = ["marker", "panel", "mash_refseq", "sourmash"] if args.mode == "all" else [args.mode]
    for mode in modes:
        try:
            if mode == "marker":
                from gside.analysis.multigene_identifier import identify_multigene

                results[mode] = identify_multigene(contigs).to_dict()
            elif mode in ("panel", "mash_refseq"):
                from gside.analysis.ani_identifier import identify_by_ani

                results[mode] = identify_by_ani(contigs, mode=mode, db_dir=args.db_dir).to_dict()
            elif mode == "sourmash":
                from gside.analysis.sourmash_identifier import identify_sourmash

                results[mode] = identify_sourmash(contigs, db_dir=args.db_dir).to_dict()
        except Exception as e:  # noqa: BLE001 — CLI 边界：单法失败不终止其他方法
            results[mode] = {"method": mode, "error": str(e)[:300]}

    if args.mode == "all" and len(results) > 1:
        verdict = _arbitrate(results)
    else:
        first = next(iter(results.values()), {})
        verdict = {
            "species": first.get("species", first.get("result", {}).get("species", "Unknown")),
            "confidence": first.get("confidence", first.get("result", {}).get("confidence", "low")),
            "basis": list(results),
        }

    payload = {
        "analysis_type": "species_identification",
        "tool": "gside",
        "version": _version(),
        "contigs": contigs,
        "methods": results,
        "verdict": verdict,
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


def _arbitrate(results: dict[str, Any]) -> dict[str, Any]:
    """ANI-layer (panel/mash/sourmash) high-confidence beats marker; else marker."""
    _LAYER = {"panel": 2, "mash_refseq": 2, "sourmash": 2, "marker": 1}
    best = None
    for mode, payload in results.items():
        if "error" in payload:
            continue
        res = payload.get("result", payload)
        species = res.get("species", "Unknown")
        conf = res.get("confidence", "low")
        score = (_LAYER.get(mode, 0), 1 if conf == "high" else 0)
        if species != "Unknown" and (best is None or score > best[0]):
            best = (score, species, conf, mode)
    if best is None:
        return {"species": "Unknown", "confidence": "low", "basis": list(results)}
    return {"species": best[1], "confidence": best[2], "basis": [best[3]]}


def _version() -> str:
    try:
        from importlib.metadata import version

        return version("gside")
    except Exception:  # noqa: BLE001
        return "0.1.0.dev"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gside", description="Genome Species IDentification Engine"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_species = sub.add_parser("species", help="Identify species from contigs")
    p_species.add_argument("contigs", help="Assembled contigs FASTA")
    p_species.add_argument(
        "--mode",
        choices=["marker", "panel", "mash_refseq", "sourmash", "all"],
        default="marker",
    )
    p_species.add_argument(
        "--db-dir", default=None, help="Database root (default: $GSIDE_DB_DIR or data/db)"
    )
    p_species.set_defaults(func=_run_species)

    parser.add_argument("--version", action="store_true")
    args = parser.parse_args(argv)
    if args.version:
        print(f"gside {_version()}")
        return 0
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
