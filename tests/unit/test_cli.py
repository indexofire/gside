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

_HAS_BIN = shutil.which("blastn") is not None or Path(
    _PROJECT, ".pixi/envs/default/bin/blastn"
).exists()

CONTIGS = _PROJECT / "data/reference/species/markers_v2.fasta"


@pytest.mark.skipif(not _HAS_BIN, reason="blastn not available")
class TestSpeciesMarker:
    def test_marker_self_hit(self):
        """markers_v2 自身作 query → 至少命中多个 marker。"""
        from gside.analysis.multigene_identifier import identify_multigene

        r = identify_multigene(str(CONTIGS)).to_dict()
        assert len(r.get("detected_markers", [])) >= 5


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
