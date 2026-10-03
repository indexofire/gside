"""Marker goldens (real blastn + tracked marker DB)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gside.config import which as _which

_REPO = Path(__file__).resolve().parents[2]
_N16961 = _REPO / "data" / "db" / "D2_ani" / "genomes" / "GCF_000006745.1_ASM674v1_genomic.fna"
_HAS_BLASTN = _which("blastn") is not None


@pytest.mark.skipif(not _HAS_BLASTN, reason="blastn not available")
@pytest.mark.skipif(not _N16961.is_file(), reason="N16961 genome missing")
class TestCholeraeN16961:
    def test_o1_toxigenic_golden(self):
        from gside.analysis.multigene_identifier import identify_multigene

        r = identify_multigene(str(_N16961))
        assert r.species == "Vibrio_cholerae_O1"
        assert r.confidence == "high"
        detected = {m["gene"] for m in r.detected_markers}
        assert {"ompw", "wben", "ctxa"} <= detected
        assert "wbfr" not in detected
        assert any("toxigenic" in note for note in r.notes)
