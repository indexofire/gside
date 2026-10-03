"""gside CLI 冒烟测试（marker 法，真实 BLAST，需 pixi 二进制）。"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PROJECT / "src"))

_BACMAP_PIXI = Path("/home/mark/repos/github/hermes-bacmap/.pixi/envs/default/bin/blastn")
if _BACMAP_PIXI.exists():
    import os

    os.environ.setdefault("GSIDE_PIXI_BIN", str(_BACMAP_PIXI.parent))

from gside.config import which as _which  # noqa: E402

_HAS_BIN = _which("blastn") is not None

CONTIGS = _PROJECT.parent / "data/db/L1_marker/markers.fasta"


@pytest.mark.skipif(not _HAS_BIN, reason="blastn not available")
class TestSpeciesMarker:
    def test_inva_single_gene_calls_salmonella(self, tmp_path):
        """提取 invA 序列作 query → 应判 Salmonella（单基因规则）。"""
        from gside.analysis.multigene_identifier import identify_multigene

        blocks = CONTIGS.read_text().split(">")
        inva = next(b for b in blocks if b.startswith("markers~~~inva~~~"))
        seq = inva.split("\n", 1)[1].replace("\n", "")
        q = tmp_path / "inva.fna"
        q.write_text(f">inva_q\n{seq}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Salmonella"

    def _gene_seq(self, prefix: str) -> str:
        blocks = CONTIGS.read_text().split(">")
        block = next(b for b in blocks if b.startswith(prefix))
        return block.split("\n", 1)[1].replace("\n", "")

    def test_excluded_gene_suppresses_coli_call(self, tmp_path):
        """ceuE + mapa → coli 被排除守卫抑制，jejuni 又凑不够 min_hits=2，应判 Unknown。"""
        from gside.analysis.multigene_identifier import identify_multigene

        ceue = self._gene_seq("markers~~~ceue~~~")
        mapa = self._gene_seq("markers~~~mapa~~~")
        q = tmp_path / "excl.fna"
        q.write_text(f">ceue_q\n{ceue}\n>mapa_q\n{mapa}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Unknown"

    def test_min_hits_two_needs_two(self, tmp_path):
        """单个 mapa 达不到 jejuni 规则 min_hits=2，应判 Unknown。"""
        from gside.analysis.multigene_identifier import identify_multigene

        mapa = self._gene_seq("markers~~~mapa~~~")
        q = tmp_path / "mapa.fna"
        q.write_text(f">mapa_q\n{mapa}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Unknown"


class TestRenderTable:
    def _payload(self):
        return {
            "contigs": "a.fna",
            "methods": {
                "marker": {
                    "species": "Salmonella",
                    "confidence": "high",
                    "matched_rule": "inva (1/1 required)",
                }
            },
            "verdict": {"species": "Salmonella", "confidence": "high", "basis": ["marker"]},
        }

    def test_tsv_single_row(self):
        from gside.cli import TABLE_COLUMNS, _render_table

        out = _render_table([self._payload()], markdown=False).splitlines()
        assert out[0].split("\t") == TABLE_COLUMNS
        row = out[1].split("\t")
        assert row[0] == "a.fna" and row[1] == "Salmonella" and row[6] == "inva (1/1 required)"

    def test_md_shape(self):
        from gside.cli import _render_table

        out = _render_table([self._payload()], markdown=True).splitlines()
        assert out[0].startswith("| contigs |") and out[1].count("---") == 14
        assert "Salmonella" in out[2]

    def test_error_cell(self):
        from gside.cli import _render_table

        payload = {
            "contigs": "b.fna",
            "methods": {"panel": {"method": "panel", "error": "db missing"}},
            "verdict": {"species": "Unknown", "confidence": "low", "basis": ["panel"]},
        }
        row = _render_table([payload], markdown=False).splitlines()[1].split("\t")
        assert row[1] == "Unknown" and "db missing" in row[13]


class TestArbitrate:
    def test_ani_beats_marker(self):
        from gside.cli import _arbitrate

        verdict = _arbitrate(
            {
                "marker": {"species": "X_marker", "confidence": "high"},
                "panel": {"result": {"species": "Y_ani", "confidence": "high"}},
            }
        )
        assert verdict["species"] == "Y_ani"
        assert verdict["basis"] == ["panel"]

    def test_errors_skipped(self):
        from gside.cli import _arbitrate

        verdict = _arbitrate(
            {"panel": {"error": "db missing"}, "marker": {"species": "M", "confidence": "high"}}
        )
        assert verdict["species"] == "M"
