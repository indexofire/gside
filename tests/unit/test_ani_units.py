"""Pure unit tests for ani_identifier (monkeypatched backends, no binaries)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gside.analysis import ani_identifier as ani
from gside.analysis.ani_identifier import (
    AniIdResult,
    _db_version,
    _hits_from_rows,
    _load_taxa_map,
    _mash_species,
    identify_by_ani,
)
from gside.engine.backends.skani import AniHit

_REPO = Path(__file__).resolve().parents[2]


@pytest.fixture()
def taxadb(tmp_path):
    (tmp_path / "metadata.tsv").write_text(
        "genome\tspecies\nGCF_A\tEscherichia coli K12\nGCF_B\tVibrio cholerae\n"
    )
    return tmp_path


class TestTaxaMap:
    def test_real_panel_map(self):
        mapping = _load_taxa_map(_REPO / "data" / "db" / "L2_ani")
        assert len(mapping) == 291
        assert mapping["GCF_000005845.2_ASM584v2_genomic.fna"].startswith("Escherichia coli")

    def test_fixture_map(self, taxadb):
        assert _load_taxa_map(taxadb)["GCF_A"] == "Escherichia coli K12"

    def test_empty_dir(self, tmp_path):
        assert _load_taxa_map(tmp_path) == {}

    def test_malformed_lines_skipped(self, tmp_path):
        (tmp_path / "metadata.tsv").write_text("genome\tspecies\nlonelyline\nA\tB\textra\n")
        assert _load_taxa_map(tmp_path) == {"A": "B"}


class TestDbVersion:
    def _db_with_sibling_manifests(self, tmp_path):
        db_dir = tmp_path / "L2_ani"
        db_dir.mkdir()
        manifests = tmp_path / "manifests"
        manifests.mkdir()
        return db_dir, manifests

    def test_missing_manifest_is_unknown(self, tmp_path):
        db_dir, _ = self._db_with_sibling_manifests(tmp_path)
        assert _db_version(db_dir, "L2_ani") == "unknown"

    def test_bad_json_is_unknown(self, tmp_path):
        db_dir, manifests = self._db_with_sibling_manifests(tmp_path)
        (manifests / "L2_ani.json").write_text("{not json")
        assert _db_version(db_dir, "L2_ani") == "unknown"

    def test_checksum_truncated(self, tmp_path):
        db_dir, manifests = self._db_with_sibling_manifests(tmp_path)
        (manifests / "L2_ani.json").write_text('{"checksum": "abcdef123456"}')
        assert _db_version(db_dir, "L2_ani") == "abcdef12"


class TestHelpers:
    def test_bad_mode(self, tmp_path):
        with pytest.raises(ValueError, match="unsupported ANI mode"):
            identify_by_ani("x.fna", mode="wat", db_dir=tmp_path)

    def test_hits_from_rows(self):
        hits = _hits_from_rows([("g1", 99.0, 0.9)])
        assert hits[0].ref == "g1" and hits[0].ani == 99.0

    def test_mash_species_passthrough(self):
        assert _mash_species("GCF_X") == "GCF_X"

    def test_to_dict_shape(self):
        d = AniIdResult(method="panel", database={"name": "L2_ani"}, result={"a": 1}).to_dict()
        assert d["analysis_type"] == "species_identification"
        assert d["result"] == {"a": 1}


def _skani_hit(ref="GCF_A", ani=99.0, af=0.9):
    return AniHit(ref=ref, ani=ani, aligned_fraction=af)


class TestPanelThresholds:
    def test_high(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_skani_search", lambda q, d: [_skani_hit(ani=99.0)])
        r = identify_by_ani("q.fna", mode="panel", db_dir=taxadb).result
        assert (r["species"], r["confidence"]) == ("Escherichia coli K12", "high")

    def test_medium_band(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_skani_search", lambda q, d: [_skani_hit(ani=94.0)])
        r = identify_by_ani("q.fna", mode="panel", db_dir=taxadb).result
        assert r["confidence"] == "medium"

    def test_low_af_rejects(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_skani_search", lambda q, d: [_skani_hit(ani=99.0, af=0.5)])
        r = identify_by_ani("q.fna", mode="panel", db_dir=taxadb).result
        assert (r["species"], r["confidence"]) == ("Unknown", "low")

    def test_below_medium(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_skani_search", lambda q, d: [_skani_hit(ani=90.0)])
        assert identify_by_ani("q.fna", mode="panel", db_dir=taxadb).result["species"] == "Unknown"

    def test_no_hits(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_skani_search", lambda q, d: [])
        r = identify_by_ani("q.fna", mode="panel", db_dir=taxadb).result
        assert r["species"] == "Unknown" and r["ani"] is None

    def test_unknown_ref_falls_back_to_id(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_skani_search", lambda q, d: [_skani_hit(ref="GCF_ZZ")])
        r = identify_by_ani("q.fna", mode="panel", db_dir=taxadb).result
        assert r["species"] == "GCF_ZZ"


class TestMashThresholds:
    def test_high(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_mash_dist", lambda q, m: [("GCF_A", 0.99)])
        r = identify_by_ani("q.fna", mode="mash_refseq", db_dir=taxadb).result
        assert (r["species"], r["confidence"]) == ("Escherichia coli K12", "high")

    def test_medium_band(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_mash_dist", lambda q, m: [("GCF_A", 0.95)])
        r = identify_by_ani("q.fna", mode="mash_refseq", db_dir=taxadb).result
        assert r["confidence"] == "medium"

    def test_low(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_mash_dist", lambda q, m: [("GCF_A", 0.5)])
        r = identify_by_ani("q.fna", mode="mash_refseq", db_dir=taxadb).result
        assert r["confidence"] == "low"

    def test_empty(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_mash_dist", lambda q, m: [])
        r = identify_by_ani("q.fna", mode="mash_refseq", db_dir=taxadb).result
        assert r["species"] == "Unknown"

    def test_bad_distance_skipped(self, monkeypatch):
        import subprocess
        import types

        def completed(stdout=""):
            return types.SimpleNamespace(returncode=0, stdout=stdout, stderr="")

        monkeypatch.setattr(subprocess, "run", lambda *a, **k: completed("refX\tq\tbadfloat"))
        assert ani._mash_dist("q.fna", "db.msh") == []
