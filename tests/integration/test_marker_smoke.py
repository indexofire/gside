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
CONTIGS = _REPO / "data" / "db" / "D1_marker" / "markers.fasta"
_HAS_BLASTN = _which("blastn") is not None


def _gene_seq(prefix: str) -> str:
    blocks = CONTIGS.read_text().split(">")
    block = next(b for b in blocks if b.startswith(prefix))
    return block.split("\n", 1)[1].replace("\n", "")


@pytest.mark.skipif(not _HAS_BLASTN, reason="blastn not available")
class TestSpeciesMarker:
    def test_inva_single_gene_calls_salmonella(self, tmp_path):
        """inva alone satisfies its single-gene rule → Salmonella_enterica."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "inva.fna"
        q.write_text(f">inva_q\n{_gene_seq('markers~~~inva~~~')}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Salmonella_enterica"

    def test_excluded_gene_suppresses_coli_call(self, tmp_path):
        """ceuE + mapa suppress the coli call; jejuni then misses min_hits=2 → Unknown."""
        from gside.analysis.multigene_identifier import identify_multigene

        ceue = _gene_seq("markers~~~ceue~~~")
        mapa = _gene_seq("markers~~~mapa~~~")
        q = tmp_path / "excl.fna"
        q.write_text(f">ceue_q\n{ceue}\n>mapa_q\n{mapa}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Unknown"

    def test_min_hits_two_needs_two(self, tmp_path):
        """A single mapa does not meet the jejuni rule's min_hits=2 → Unknown."""
        from gside.analysis.multigene_identifier import identify_multigene

        q = tmp_path / "mapa.fna"
        q.write_text(f">mapa_q\n{_gene_seq('markers~~~mapa~~~')}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Unknown"
