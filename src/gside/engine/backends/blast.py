from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from .._env import require_bin, run_checked, which
from ..hits import Hit

_BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore qlen slen"
)

_PARAM_MAP = {
    "pident": "perc_identity",
    "perc_identity": "perc_identity",
    "qcov_hsp_perc": "qcov_hsp_perc",
    "threads": "num_threads",
    "num_threads": "num_threads",
    "evalue": "evalue",
    "max_targets": "max_target_seqs",
    "max_target_seqs": "max_target_seqs",
}


class BlastBackend:
    """BLAST+ backend supporting blastn, blastp, blastx, tblastn, tblastx."""

    def __init__(self, tool: str = "blastn", threads: int = 4, binary: str | None = None):
        self.tool = tool
        self.threads = threads
        # binary= (caller-resolved) skips discovery so callers keep their own error messages.
        self._bin = binary if binary is not None else self._find_binary()

    def _find_binary(self) -> str:
        return require_bin(self.tool, hint=f"{self.tool} not found in PATH", resolver=which)

    def make_db(self, fasta_file: Path, db_path: Path, db_type: str = "nucl") -> None:
        makeblastdb = which("makeblastdb")
        if not makeblastdb:
            raise RuntimeError("makeblastdb not found")
        cmd = [makeblastdb, "-in", str(fasta_file), "-dbtype", db_type, "-out", str(db_path)]
        subprocess.run(cmd, check=True, capture_output=True, timeout=120)

    def ensure_index(self, db_prefix: str, db_type: str = "nucl") -> None:
        ext = ".phr" if db_type == "prot" else ".nhr"
        if not Path(f"{db_prefix}{ext}").exists():
            fasta_candidates = [
                Path(f"{db_prefix}.fasta"),
                Path(f"{db_prefix}.fna"),
                Path(f"{db_prefix}.fa"),
                Path(f"{db_prefix}_sequences.fasta"),
                Path(f"{db_prefix}_abricate.fasta"),
            ]
            for fc in fasta_candidates:
                if fc.exists():
                    self.make_db(fc, Path(db_prefix), db_type)
                    return
            raise FileNotFoundError(f"No BLAST index or source FASTA for {db_prefix}")

    def find(
        self,
        query: Path,
        db_path: str,
        min_identity: float = 0.0,
        min_coverage: float = 0.0,
        evalue: float = 1e-5,
        max_targets: int = 500,
        **kwargs: Any,
    ) -> list[Hit]:
        cmd = [
            self._bin,
            "-query",
            str(query),
            "-db",
            db_path,
            "-outfmt",
            _BLAST_OUTFMT,
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            str(max_targets),
            "-num_threads",
            str(self.threads),
        ]

        for key, value in kwargs.items():
            if value is None or key in ("threads",):
                continue
            if key == "num_threads":
                cmd[cmd.index("-num_threads") + 1] = str(value)
                continue
            mapped = _PARAM_MAP.get(key, key)
            if isinstance(value, bool):
                if value:
                    cmd.append(f"-{mapped}")
            else:
                cmd.extend([f"-{mapped}", str(value)])

        proc = run_checked(cmd, timeout=600, name=f"{self.tool} failed")

        hits: list[Hit] = []
        for line in proc.stdout.splitlines():
            if not line.strip():
                continue
            try:
                hit = Hit.from_blast_line(line)
            except ValueError:
                continue
            if hit.identity >= min_identity and hit.query_coverage >= min_coverage:
                hits.append(hit)

        return hits
