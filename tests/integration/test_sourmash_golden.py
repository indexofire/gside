"""Sourmash binary goldens (real provisioned binary, tiny inputs, no GTDB db)."""

from __future__ import annotations

import shutil

import pytest

_has_sourmash = shutil.which("sourmash") is not None


@pytest.mark.skipif(not _has_sourmash, reason="sourmash not available")
class TestSourmashSketch:
    def test_sketch_and_run(self, tmp_path):
        import random

        from gside.analysis.sourmash_identifier import _run

        random.seed(42)
        fasta = tmp_path / "tiny.fna"
        fasta.write_text(">s1\n" + "".join(random.choice("ACGT") for _ in range(20000)) + "\n")
        sig = tmp_path / "tiny.sig"
        _run(["sourmash", "sketch", "dna", "-p", "k=31,scaled=1000", str(fasta), "-o", str(sig)])
        assert sig.is_file()

    def test_failure_raises(self):
        from gside.analysis.sourmash_identifier import _run

        with pytest.raises(RuntimeError, match="failed"):
            _run(["sourmash", "gather", "nonexistent.sig", "nonexistent.zip"])
