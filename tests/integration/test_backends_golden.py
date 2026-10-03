"""Backend integration goldens (real provisioned binaries, tiny inputs)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]

_has_mmseqs = shutil.which("mmseqs") is not None


@pytest.mark.skipif(not _has_mmseqs, reason="mmseqs not available")
class TestMmseqsLinclust:
    def test_identical_pair_clusters(self, tmp_path):
        from gside.engine.backends.mmseqs2 import Mmseqs2Backend, parse_cluster_tsv

        fasta = tmp_path / "seqs.fna"
        fasta.write_text(
            ">seq1\n" + "ACGT" * 50 + "\n>seq2\n" + "ACGT" * 50 + "\n>seq3\n" + "TGCA" * 50 + "\n"
        )
        tsv = Mmseqs2Backend().cluster(fasta, tmp_path / "out")
        clusters = parse_cluster_tsv(tsv)
        pair = next(v for v in clusters.values() if "seq1" in v)
        assert "seq2" in pair


def _has_kma() -> bool:
    from gside.config import which

    return which("kma") is not None


class TestKmaMapReads:
    def test_index_and_map(self, tmp_path):
        import pytest

        if not _has_kma():
            pytest.skip("kma not available")
        from gside.engine.backends.kma import KmaBackend

        tpl = tmp_path / "tpl.fna"
        tpl.write_text(">tpl1\n" + "ACGT" * 125 + "\n")
        reads = tmp_path / "r1.fq"
        read = "ACGT" * 38
        reads.write_text(f"@r1\n{read}\n+\n{'I' * len(read)}\n")
        backend = KmaBackend()
        backend.make_index(tpl, tmp_path / "idx")
        hits = backend.find(reads, tmp_path / "idx", output_dir=tmp_path / "kma_out")
        assert hits and hits[0].subject_id == "tpl1"
        assert hits[0].backend == "kma"
