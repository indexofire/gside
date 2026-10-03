from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .._env import require_bin, run_checked, which


@dataclass(frozen=True)
class AniHit:
    ref: str
    ani: float
    aligned_fraction: float
    backend: str = "skani"


def _find_bin() -> str:
    return require_bin(
        "skani",
        hint=(
            "skani not found in PATH. Install: pixi add 'skani=0.3.*' (locked to "
            "match the official pre-sketched GTDB database format)"
        ),
        resolver=which,
    )


class SkaniBackend:
    """skani ANI search against a pre-sketched database (skani sketch ... -o db)."""

    def __init__(self) -> None:
        self._bin = _find_bin()

    def search(self, query: Path, db: Path) -> list[AniHit]:
        result = run_checked(
            [self._bin, "search", str(query), "-d", str(db)],
            timeout=600,
            name="skani search failed",
            exit_in_msg=False,
            truncate=None,
        )

        lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
        if not lines:
            return []

        header = [h.strip() for h in lines[0].split("\t")]
        ref_i = next((i for i, h in enumerate(header) if h in ("ref_filename", "Ref_file")), None)
        qaf_i = next(
            (
                i
                for i, h in enumerate(header)
                if h in ("Estimated_query_aligned_fraction", "Align_fraction_query")
            ),
            None,
        )
        ani_i = next((i for i, h in enumerate(header) if h == "ANI"), None)
        if ref_i is None or qaf_i is None or ani_i is None:
            ref_i, qaf_i, ani_i = 0, 3, 5

        hits: list[AniHit] = []
        for ln in lines[1:]:
            cols = ln.split("\t")
            if len(cols) <= max(ref_i, qaf_i, ani_i):
                continue
            hits.append(
                AniHit(
                    ref=cols[ref_i].strip(),
                    ani=float(cols[ani_i]),
                    aligned_fraction=float(cols[qaf_i]),
                )
            )
        return hits
