"""gside CLI 冒烟测试（marker 法，真实 BLAST，需 pixi 二进制）。"""

from __future__ import annotations

import json
import shutil
import subprocess
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

CONTIGS = _PROJECT.parent / "data/reference/species/markers_v2.fasta"


@pytest.mark.skipif(not _HAS_BIN, reason="blastn not available")
class TestSpeciesMarker:
    def test_inva_single_gene_calls_salmonella(self, tmp_path):
        """提取 invA 序列作 query → 应判 Salmonella（单基因规则）。"""
        from gside.analysis.multigene_identifier import identify_multigene

        blocks = CONTIGS.read_text().split(">")
        inva = next(b for b in blocks if b.startswith("markers_v2~~~inva~~~"))
        seq = inva.split("\n", 1)[1].replace("\n", "")
        q = tmp_path / "inva.fna"
        q.write_text(f">inva_q\n{seq}\n")
        r = identify_multigene(str(q)).to_dict()
        assert r["species"] == "Salmonella"


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
