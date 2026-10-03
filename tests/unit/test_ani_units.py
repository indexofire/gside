"""Pure unit tests for ani_identifier (monkeypatched backends, no binaries)."""

from __future__ import annotations

import subprocess
import types
from pathlib import Path

import pytest

from gside import config
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


def _completed(stdout: str = "", returncode: int = 0, stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


@pytest.fixture()
def taxadb(tmp_path):
    (tmp_path / "metadata.tsv").write_text(
        "genome\tspecies\nGCF_A\tEscherichia coli K12\nGCF_B\tVibrio cholerae\n"
    )
    (tmp_path / "mash.msh").write_bytes(b"")
    return tmp_path


class TestTaxaMap:
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
        monkeypatch.setattr(config, "which", lambda tool: "/fake/mash")
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("refX\tq\tbadfloat"))
        assert ani._mash_dist("q.fna", "db.msh") == []


class TestMashDelegation:
    def test_routes_through_mash_backend(self, monkeypatch, tmp_path):
        from gside.engine.backends import kmer as kmer_mod

        calls: dict[str, object] = {}

        class _KD:
            def __init__(self, ref: str, dist: float) -> None:
                self.reference_id = ref
                self.distance = dist

        class _M:
            def __init__(self, kmer_size=21, sketch_size=1000, threads=4, binary=None):
                calls["ctor"] = binary

            def distance(self, query, reference, max_distance=0.1, **kw):
                calls["find"] = (query, reference, max_distance)
                return [_KD("GCF_A", 0.5), _KD("GCF_B", 0.25)]

        monkeypatch.setattr(kmer_mod, "MashBackend", _M)
        monkeypatch.setattr(config, "which", lambda tool: "/fake/mash")
        rows = ani._mash_dist("q.fna", tmp_path / "mash.msh")
        assert rows == [("GCF_A", 0.5), ("GCF_B", 0.75)]
        assert calls["ctor"] == "/fake/mash"
        assert calls["find"] == (tmp_path / "mash.msh", Path("q.fna"), 1.0)


class TestMashFailures:
    def test_missing_db_raises_with_setup_hint(self, tmp_path, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(""))
        with pytest.raises(RuntimeError, match="gside db setup --tier mash"):
            identify_by_ani("q.fna", mode="mash_refseq", db_dir=tmp_path)

    def test_missing_binary_raises_with_fix_hint(self, monkeypatch, tmp_path):
        (tmp_path / "mash.msh").write_bytes(b"")
        monkeypatch.setattr(config, "which", lambda tool: None)
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(""))
        with pytest.raises(RuntimeError, match=r"mash not found.*GSIDE_PIXI_BIN"):
            ani._mash_dist("q.fna", tmp_path / "mash.msh")

    def test_nonzero_exit_raises_runtime_error(self, monkeypatch, tmp_path):
        (tmp_path / "mash.msh").write_bytes(b"")
        monkeypatch.setattr(config, "which", lambda tool: "/fake/mash")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: _completed("", returncode=1, stderr="mash exploded"),
        )
        with pytest.raises(RuntimeError, match=r"mash dist failed \(exit 1\)"):
            ani._mash_dist("q.fna", tmp_path / "mash.msh")

    def test_failure_message_carries_stderr(self, monkeypatch, tmp_path):
        (tmp_path / "mash.msh").write_bytes(b"")
        monkeypatch.setattr(config, "which", lambda tool: "/fake/mash")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: _completed("", returncode=2, stderr="sketch version mismatch"),
        )
        with pytest.raises(RuntimeError, match="sketch version mismatch"):
            ani._mash_dist("q.fna", tmp_path / "mash.msh")
