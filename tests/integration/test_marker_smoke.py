"""Marker-mode smoke tests: real blastn against the tracked marker database.

These are integration tests on purpose — they need the conda-provisioned
``blastn`` binary (resolved via ``GSIDE_PIXI_BIN`` or the pixi env). The
hermetic unit suite lives in ``tests/unit/``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gside.config import which as _which

_REPO = Path(__file__).resolve().parents[2]
CONTIGS = _REPO / "data" / "db" / "D1_marker" / "default" / "markers.fasta"
_HAS_BLASTN = _which("blastn") is not None


def _gene_seq(prefix: str) -> str:
    blocks = CONTIGS.read_text().split(">")
    block = next(b for b in blocks if b.startswith(prefix))
    return block.split("\n", 1)[1].replace("\n", "")


@pytest.mark.skipif(not _HAS_BLASTN, reason="blastn not available")
class TestSpeciesMarker:
    def test_inva_alone_falls_back_to_genus(self, tmp_path):
        """inva alone fires only the genus fallback (species needs the anchor + marker)."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "inva.fna"
        q.write_text(f">inva_q\n{_gene_seq('gside|db=species_markers_v2|gene=inva|acc=')}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Salmonella"

    def test_inva_plus_spi2_calls_enterica(self, tmp_path):
        """inva (required anchor) + ssav (SPI-2 marker) → Salmonella_enterica."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "inva_ssav.fna"
        q.write_text(
            f">inva_q\n{_gene_seq('gside|db=species_markers_v2|gene=inva|acc=')}\n"
            f">ssav_q\n{_gene_seq('gside|db=species_markers_v2|gene=ssav|acc=')}\n"
        )
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Salmonella_enterica"

    def test_spi2_without_inva_calls_nothing(self, tmp_path):
        """SPI-2 markers without the inva anchor never fire the species rule."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "spi2.fna"
        q.write_text(
            f">ssav_q\n{_gene_seq('gside|db=species_markers_v2|gene=ssav|acc=')}\n"
            f">ssec_q\n{_gene_seq('gside|db=species_markers_v2|gene=ssec|acc=')}\n"
        )
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Unknown"

    def test_chromosomal_anchor_calls_anthracis(self, tmp_path):
        """adk (group anchor) + ba5345 (lambdaBa03) → Bacillus_anthracis —
        plasmid-cured/vaccine strains (Sterne-like) stay anthracis, and the
        plasmid genes never gate the species call. Regression guard for the
        corrupted-marker bug (the original paga/capb entries matched no
        B. anthracis genome, so the rule could never fire)."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "ba.fna"
        q.write_text(
            f">adk_q\n{_gene_seq('gside|db=species_markers_v2|gene=adk|acc=')}\n"
            f">ba5345_q\n{_gene_seq('gside|db=species_markers_v2|gene=ba5345|acc=')}\n"
            f">paga_q\n{_gene_seq('gside|db=species_markers_v2|gene=paga|acc=')}\n"
            f">capb_q\n{_gene_seq('gside|db=species_markers_v2|gene=capb|acc=')}\n"
        )
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Bacillus_anthracis"
        assert {"adk", "ba5345", "paga", "capb"} <= {m["gene"] for m in r["detected_markers"]}

        # 无质毒基因（疫苗株画像）仍是炭疽
        q2 = tmp_path / "ba_sterne.fna"
        q2.write_text(
            f">adk_q\n{_gene_seq('gside|db=species_markers_v2|gene=adk|acc=')}\n"
            f">ba5345_q\n{_gene_seq('gside|db=species_markers_v2|gene=ba5345|acc=')}\n"
        )
        r2 = identify_multigene(str(q2)).to_dict()
        assert r2["species"] == "Bacillus_anthracis"

    def test_plasmids_without_chromosomal_anchor_stay_group(self, tmp_path):
        """Bcbva/cross-over profile: adk + paga/capb, no ba5345 →
        Bacillus_cereus_group with the anthrax-like note."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "bcbva.fna"
        q.write_text(
            f">adk_q\n{_gene_seq('gside|db=species_markers_v2|gene=adk|acc=')}\n"
            f">paga_q\n{_gene_seq('gside|db=species_markers_v2|gene=paga|acc=')}\n"
            f">capb_q\n{_gene_seq('gside|db=species_markers_v2|gene=capb|acc=')}\n"
        )
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Bacillus_cereus_group"
        assert any("Anthrax-like" in n for n in r["notes"])

    def test_adk_alone_calls_cereus_group(self, tmp_path):
        """Plain cereus/thuringiensis (adk only) → Bacillus_cereus_group."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "cereus.fna"
        q.write_text(f">adk_q\n{_gene_seq('gside|db=species_markers_v2|gene=adk|acc=')}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Bacillus_cereus_group"

    def test_excluded_gene_suppresses_coli_call(self, tmp_path):
        """ceuE + hipo suppress the coli call; jejuni then misses min_hits=2 → Unknown."""
        from gside.analysis.multigene_identifier import identify_multigene

        ceue = _gene_seq("gside|db=species_markers_v2|gene=ceue|acc=")
        hipo = _gene_seq("gside|db=species_markers_v2|gene=hipo|acc=")
        q = tmp_path / "excl.fna"
        q.write_text(f">ceue_q\n{ceue}\n>hipo_q\n{hipo}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Unknown"

    def test_mapa_no_longer_excludes_coli(self, tmp_path):
        """ceuE + mapa fires the coli call (mapa cross-reacts with C. coli)."""
        from gside.analysis.multigene_identifier import identify_multigene

        ceue = _gene_seq("gside|db=species_markers_v2|gene=ceue|acc=")
        mapa = _gene_seq("gside|db=species_markers_v2|gene=mapa|acc=")
        q = tmp_path / "mapa.fna"
        q.write_text(f">ceue_q\n{ceue}\n>mapa_q\n{mapa}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Campylobacter_coli"

    def test_min_hits_two_needs_two(self, tmp_path):
        """A single mapa does not meet the jejuni rule's min_hits=2 → Unknown."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "mapa.fna"
        q.write_text(f">mapa_q\n{_gene_seq('gside|db=species_markers_v2|gene=mapa|acc=')}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Unknown"

    def test_chromosomal_anchor_calls_pseudomallei(self, tmp_path):
        """glgB anchor calls B. pseudomallei; BimA/T3SS are reporters only.

        Regression guard for the corrupted tts1/bimabp entries, which matched
        0/100 clinical B. pseudomallei genomes. The anchor must fire alone
        (BimA is absent from ~5% of genomes, one panel genome lost the TTS
        cluster) and the reporters must not fire without it."""
        from gside.analysis.multigene_identifier import identify_multigene

        glgb = _gene_seq("gside|db=species_markers_v2|gene=glgb-bp|acc=")
        bima = _gene_seq("gside|db=species_markers_v2|gene=bima-bp|acc=")
        sctn = _gene_seq("gside|db=species_markers_v2|gene=sctn-bp|acc=")

        anchor = tmp_path / "bp_anchor.fna"
        anchor.write_text(f">glgb_q\n{glgb}\n")
        r = identify_multigene(str(anchor)).to_dict()
        assert r["species"] == "Burkholderia_pseudomallei"

        full = tmp_path / "bp_full.fna"
        full.write_text(f">glgb_q\n{glgb}\n>bima_q\n{bima}\n>sctn_q\n{sctn}\n")
        r = identify_multigene(str(full)).to_dict()
        assert r["species"] == "Burkholderia_pseudomallei"
        assert {"glgb-bp", "bima-bp", "sctn-bp"} <= {m["gene"] for m in r["detected_markers"]}

        reporters = tmp_path / "bp_reporters.fna"
        reporters.write_text(f">bima_q\n{bima}\n>sctn_q\n{sctn}\n")
        r = identify_multigene(str(reporters)).to_dict()
        assert r["species"] == "Unknown"

    def test_shigella_genus_call_points_to_gspec(self, tmp_path):
        """uidA+ipaH (lacY-negative) reports genus Shigella with a gspec note."""
        from gside.analysis.multigene_identifier import identify_multigene

        uida = _gene_seq("gside|db=species_markers_v2|gene=uida|acc=")
        ipah = _gene_seq("gside|db=species_markers_v2|gene=ipah|acc=")
        q = tmp_path / "shigella.fna"
        q.write_text(f">uida_q\n{uida}\n>ipah_q\n{ipah}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Shigella"
        assert r["confidence"] == "high"
        assert any("gspec" in note.lower() for note in r["notes"])
