"""Panel/mash self-identity goldens (real binaries + real databases)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
_PANEL = _REPO / "data" / "db" / "D2_ani"
_MASH = _REPO / "data" / "db" / "D3_mash"
_K12 = _PANEL / "genomes" / "GCF_000005845.2_ASM584v2_genomic.fna"

_has_skani = shutil.which("skani") is not None
_has_mash = shutil.which("mash") is not None


@pytest.mark.skipif(not _K12.is_file(), reason="panel genome missing")
@pytest.mark.skipif(not (_PANEL / "panel.sketch").exists(), reason="panel sketch missing")
@pytest.mark.skipif(not _has_skani, reason="skani not available")
class TestPanelSelfIdentity:
    def test_k12_against_own_panel(self):
        from gside.analysis.ani_identifier import identify_by_ani

        r = identify_by_ani(str(_K12), mode="panel").result
        assert "Escherichia coli" in r["species"]
        assert r["confidence"] == "high"
        assert r["ani"] >= 99.0
        assert r["aligned_fraction"] >= 0.65
        assert r["top_hits"][0]["genome"] == _K12.name


@pytest.mark.skipif(not _K12.is_file(), reason="panel genome missing")
@pytest.mark.skipif(
    not (_MASH / "mash.msh").exists() and not (_MASH / "payload.bin").exists(),
    reason="mash sketch missing",
)
@pytest.mark.skipif(not _has_mash, reason="mash not available")
class TestMashSelfIdentity:
    def test_k12_against_own_sketch(self):
        from gside.analysis.ani_identifier import identify_by_ani

        r = identify_by_ani(str(_K12), mode="mash_refseq").result
        assert r["top_hits"], "expected at least one mash hit"
        assert 0.0 <= r["identity"] <= 1.0


@pytest.mark.skipif(not (_PANEL / "metadata.tsv").is_file(), reason="panel metadata missing")
class TestPanelTaxaMap:
    def test_real_panel_map_has_291_genomes(self):
        from gside.analysis.ani_identifier import _load_taxa_map

        mapping = _load_taxa_map(_PANEL)
        assert len(mapping) == 291
        assert mapping["GCF_000005845.2_ASM584v2_genomic.fna"].startswith("Escherichia coli")


_MASH_GOLDENS = {
    "GCF_000006945.2_ASM694v2_genomic.fna": ("Salmonella_enterica", 1.0),
    "GCF_000009085.1_ASM908v1_genomic.fna": ("Campylobacter_jejuni", 1.0),
    "GCF_000465235.1_ASM46523v1_genomic.fna": ("Campylobacter_coli", 0.99),
    "GCF_000005845.2_ASM584v2_genomic.fna": ("Escherichia_coli", 1.0),
    "GCF_000009005.1_ASM900v1_genomic.fna": ("Staphylococcus_aureus", 0.98),
}


@pytest.mark.skipif(not _has_mash, reason="mash not available")
class TestMashGoldens:
    @pytest.mark.parametrize("fna,expected", list(_MASH_GOLDENS.items()))
    def test_mash_golden(self, fna, expected):
        from gside.analysis.ani_identifier import identify_by_ani

        genome = _REPO / "data" / "db" / "D2_ani" / "genomes" / fna
        if not genome.is_file():
            pytest.skip("genome missing")
        want_species, want_identity = expected
        r = identify_by_ani(str(genome), mode="mash_refseq").result
        assert r["species"].startswith(want_species)
        assert r["confidence"] == "high"
        assert r["identity"] >= want_identity
