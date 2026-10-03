"""Marker-vs-panel agreement goldens on panel's own genomes (real binaries)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
_GENOMES = _REPO / "data" / "db" / "L2_ani" / "genomes"
_SKETCH = _REPO / "data" / "db" / "L2_ani" / "panel.sketch"

_CASES = {
    "GCF_000006945.2_ASM694v2_genomic.fna": ("Salmonella", "Salmonella"),
    "GCF_000009085.1_ASM908v1_genomic.fna": ("Campylobacter_jejuni", "Campylobacter jejuni"),
    "GCF_000465235.1_ASM46523v1_genomic.fna": ("Campylobacter_coli", "Campylobacter coli"),
    "GCF_000005845.2_ASM584v2_genomic.fna": ("DEC", "Escherichia coli"),
}

_ready = (
    _SKETCH.exists()
    and all((_GENOMES / f).is_file() for f in _CASES)
    and shutil.which("blastn") is not None
    and shutil.which("skani") is not None
)


@pytest.mark.skipif(not _ready, reason="panel genomes/binaries missing")
class TestMarkerPanelAgreement:
    @pytest.mark.parametrize("fna,expected", list(_CASES.items()))
    def test_agreement(self, fna, expected):
        from gside.analysis.ani_identifier import identify_by_ani
        from gside.analysis.multigene_identifier import identify_multigene

        marker_species, panel_fragment = expected
        m = identify_multigene(str(_GENOMES / fna)).to_dict()
        p = identify_by_ani(str(_GENOMES / fna), mode="panel").result
        assert m["species"] == marker_species
        assert m["confidence"] == "high"
        assert panel_fragment in p["species"]
        assert p["confidence"] == "high"
