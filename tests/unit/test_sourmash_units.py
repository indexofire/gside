"""Pure unit tests for sourmash_identifier (backends monkeypatched)."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from gside.analysis import sourmash_identifier as sm
from gside.analysis.sourmash_identifier import (
    SourmashIdResult,
    _db_version,
    _lineage_map,
    identify_by_sourmash,
)

_GATHER_CSV = (
    "intersect_bp,f_match,f_unique_weighted,filename,name\n900,1.0,0.9,/tmp/db.zip,GCF_A\n"
)


def _seed_db(db_dir: Path) -> None:
    (db_dir / "gtdb-reps-k31.zip").write_bytes(b"x")
    (db_dir / "lineages.csv").write_text("GCF_A,s__Vp\n")


def _fake_run_writing_outputs(csv_text: str = _GATHER_CSV):
    """Fake _run mimicking sourmash: writes whatever path follows ``-o``."""

    def run(cmd: list[str]) -> str:
        out = Path(cmd[cmd.index("-o") + 1])
        if cmd[1] == "sketch":
            out.write_bytes(b"sig")
        else:
            out.write_text(csv_text)
        return ""

    return run


class TestLineageMap:
    def test_last_rank_wins(self, tmp_path):
        lin = tmp_path / "lineages.csv"
        lin.write_text(
            "GCF_A,d__Bacteria;p__Proteobacteria;c__Gamma;s__Vibrio parahaemolyticus\nshort\n"
        )
        assert _lineage_map(lin) == {"GCF_A": "s__Vibrio parahaemolyticus"}

    def test_missing_file(self, tmp_path):
        assert _lineage_map(tmp_path / "nope.csv") == {}


class TestDbVersion:
    def test_default_rs226(self, tmp_path):
        assert _db_version(tmp_path / "L4_sourmash") == "RS226"

    def test_checksum(self, tmp_path):
        (tmp_path / "L4_sourmash").mkdir()
        manifests = tmp_path / "manifests"
        manifests.mkdir()
        (manifests / "L4_sourmash.json").write_text('{"checksum": "abcdef123456"}')
        assert _db_version(tmp_path / "L4_sourmash") == "abcdef12"


class TestToDict:
    def test_shape(self):
        d = SourmashIdResult(
            method="sourmash", database={"name": "L4_sourmash"}, result={"a": 1}
        ).to_dict()
        assert d["analysis_type"] == "species_identification"
        assert d["result"] == {"a": 1}


def _partition(*members):
    return [{"lineage": lin, "f_unique_weighted": fuw} for lin, fuw in members]


class TestConfidence:
    def test_high(self, monkeypatch, tmp_path):
        monkeypatch.setattr(sm, "_gather_and_tax", lambda c, b: _partition(("Vp", 0.95)))
        r = identify_by_sourmash("q.fna", db_dir=tmp_path).result
        assert (r["species"], r["confidence"]) == ("Vp", "high")
        assert r["flags"] == []

    def test_mixture_flag_blocks_high(self, monkeypatch, tmp_path):
        monkeypatch.setattr(
            sm,
            "_gather_and_tax",
            lambda c, b: _partition(("Vp", 0.95), ("Vc", 0.20)),
        )
        r = identify_by_sourmash("q.fna", db_dir=tmp_path).result
        assert r["flags"] == ["possible_mixture"]
        assert r["confidence"] == "medium"

    def test_medium(self, monkeypatch, tmp_path):
        monkeypatch.setattr(sm, "_gather_and_tax", lambda c, b: _partition(("Vp", 0.8)))
        assert identify_by_sourmash("q.fna", db_dir=tmp_path).result["confidence"] == "medium"

    def test_low_empty(self, monkeypatch, tmp_path):
        monkeypatch.setattr(sm, "_gather_and_tax", lambda c, b: [])
        r = identify_by_sourmash("q.fna", db_dir=tmp_path).result
        assert (r["species"], r["confidence"]) == ("Mixed/Unknown", "low")


class TestGatherAndTax:
    def test_csv_parsing(self, monkeypatch, tmp_path):
        import gside.analysis.sourmash_identifier as mod

        _seed_db(tmp_path)
        csv_text = (
            "intersect_bp,f_match,f_unique_weighted,filename,name\n900,1.0,0.9,/tmp/db.zip,GCF_A\n"
        )
        monkeypatch.setattr(mod, "_run", _fake_run_writing_outputs(csv_text))
        out = mod._gather_and_tax("q.fna", tmp_path)
        assert out[0] == {"lineage": "s__Vp", "f_unique_weighted": 0.9}

    def test_missing_db_files(self, tmp_path):
        import gside.analysis.sourmash_identifier as mod

        with pytest.raises(RuntimeError, match="not found"):
            mod._gather_and_tax("q.fna", tmp_path)

    def test_gather_no_match_is_empty(self, monkeypatch, tmp_path):
        import gside.analysis.sourmash_identifier as mod

        _seed_db(tmp_path)

        def fake_run(cmd):
            if cmd[1] == "gather":
                raise RuntimeError("sourmash gather failed: no matches")

        monkeypatch.setattr(mod, "_run", fake_run)
        assert mod._gather_and_tax("q.fna", tmp_path) == []

    def test_bad_fuw_is_zero(self, monkeypatch, tmp_path):
        import gside.analysis.sourmash_identifier as mod

        _seed_db(tmp_path)
        bad_csv = (
            "intersect_bp,f_match,f_unique_weighted,filename,name\n"
            "900,1.0,notanum,/tmp/db.zip,GCF_A\n"
        )
        monkeypatch.setattr(mod, "_run", _fake_run_writing_outputs(bad_csv))
        out = mod._gather_and_tax("q.fna", tmp_path)
        assert out[0] == {"lineage": "s__Vp", "f_unique_weighted": 0.0}

    def test_bad_json_version(self, tmp_path):
        (tmp_path / "L4_sourmash").mkdir()
        manifests = tmp_path / "manifests"
        manifests.mkdir()
        (manifests / "L4_sourmash.json").write_text("{not json")
        assert _db_version(tmp_path / "L4_sourmash") == "RS226"


class TestMissingBinary:
    def test_binary_absence_propagates(self, monkeypatch, tmp_path):
        (tmp_path / "gtdb-reps-k31.zip").write_bytes(b"x")
        (tmp_path / "lineages.csv").write_text("GCF_A,s__Vp\n")
        monkeypatch.setattr(
            sm, "_run", lambda cmd: (_ for _ in ()).throw(FileNotFoundError("sourmash"))
        )
        with pytest.raises(FileNotFoundError):
            identify_by_sourmash("q.fna", db_dir=tmp_path)


class TestTempDirIntermediates:
    def test_no_intermediates_left_in_db_dir(self, monkeypatch, tmp_path):
        _seed_db(tmp_path)
        outs: list[Path] = []

        def fake_run(cmd):
            out = Path(cmd[cmd.index("-o") + 1])
            outs.append(out)
            if cmd[1] == "sketch":
                out.write_bytes(b"sig")
            else:
                out.write_text(_GATHER_CSV)
            return ""

        monkeypatch.setattr(sm, "_run", fake_run)
        result = identify_by_sourmash("q.fna", db_dir=tmp_path)
        assert result.result["species"] == "s__Vp"
        assert not (tmp_path / "query.sig").exists()
        assert not (tmp_path / "gather.csv").exists()
        assert outs, "expected sketch and gather to run"
        assert all(not out.is_relative_to(tmp_path) for out in outs)

    def test_temp_dir_created_and_cleaned_up(self, monkeypatch, tmp_path):
        _seed_db(tmp_path)
        real_tmpdir = tempfile.TemporaryDirectory
        created: list[tempfile.TemporaryDirectory] = []

        def recording_tmpdir(*args, **kwargs):
            handle = real_tmpdir(*args, **kwargs)
            created.append(handle)
            return handle

        monkeypatch.setattr(tempfile, "TemporaryDirectory", recording_tmpdir)
        monkeypatch.setattr(sm, "_run", _fake_run_writing_outputs())

        identify_by_sourmash("q.fna", db_dir=tmp_path)

        assert created, "intermediates must go through tempfile.TemporaryDirectory"
        for handle in created:
            assert not Path(handle.name).exists(), f"temp dir not cleaned up: {handle.name}"
