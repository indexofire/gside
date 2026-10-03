from __future__ import annotations

from pathlib import Path

from .._env import require_bin, run_checked, which


def _find_bin() -> str:
    return require_bin(
        "mmseqs",
        hint=(
            "mmseqs not found in PATH. Install: pixi add 'mmseqs2>=15' "
            "(linear-time clustering engine for pan-genome discovery)"
        ),
        resolver=which,
    )


def parse_cluster_tsv(tsv: Path) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = {}
    for line in tsv.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 2:
            continue
        rep, member = cols[0].strip(), cols[1].strip()
        if not rep or not member:
            continue
        clusters.setdefault(rep, []).append(member)
    return clusters


class Mmseqs2Backend:
    """mmseqs2 sequence clustering via easy-linclust (linear-time)."""

    def __init__(self, threads: int = 4) -> None:
        self._bin = _find_bin()
        self._threads = threads

    def cluster(
        self,
        fasta: Path,
        out_prefix: Path,
        min_seq_id: float = 0.9,
        coverage: float = 0.8,
    ) -> Path:
        result = run_checked(
            [
                self._bin,
                "easy-linclust",
                str(fasta),
                str(out_prefix),
                str(out_prefix.parent / "mmseqs_tmp"),
                "--min-seq-id",
                str(min_seq_id),
                "-c",
                str(coverage),
                "--threads",
                str(self._threads),
            ],
            timeout=3600,
            name="mmseqs easy-linclust failed",
            exit_in_msg=False,
        )
        tsv = Path(f"{out_prefix}_cluster.tsv")
        if not tsv.exists():
            raise RuntimeError(
                f"mmseqs easy-linclust did not produce {tsv} "
                f"(stdout tail: {result.stdout.strip()[-300:]})"
            )
        return tsv
