"""Pure unit tests for read_mapper helpers (no binaries)."""

from __future__ import annotations

import subprocess
import types

import pytest

from gside.engine import read_mapper as rm
from gside.engine.read_mapper import (
    BwaReadMapper,
    Minimap2ReadMapper,
    ReadMapper,
    _get_mapper,
    _sniff_read_type,
)


def _completed(stdout="", returncode=0, stderr=""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


class TestEnsureBwaIndex:
    def test_existing_index_noop(self, tmp_path, monkeypatch):
        ref = tmp_path / "ref.fna"
        ref.write_text(">a\nACGT\n")
        (tmp_path / "ref.fna.bwt").write_bytes(b"x")
        monkeypatch.setattr(rm, "which", lambda tool: None)
        assert rm._ensure_bwa_index(str(ref)) is None

    def test_missing_binary(self, tmp_path, monkeypatch):
        monkeypatch.setattr(rm, "which", lambda tool: None)
        with pytest.raises(RuntimeError, match="bwa not found"):
            rm._ensure_bwa_index(str(tmp_path / "ref.fna"))


class TestGetMapper:
    def test_names(self):
        assert isinstance(_get_mapper("bwa"), BwaReadMapper)
        assert isinstance(_get_mapper(" BWA-MEM "), BwaReadMapper)
        assert isinstance(_get_mapper("Minimap2"), Minimap2ReadMapper)

    def test_unknown(self):
        with pytest.raises(KeyError):
            _get_mapper("nope")


class TestMapValidation:
    def test_bad_read_type(self, tmp_path):
        with pytest.raises(ValueError, match="read_type must be"):
            ReadMapper.map(["a.fq"], "ref", "out.bam", read_type="wat")


class TestAlignPaths:
    def test_no_samtools(self, monkeypatch, tmp_path):
        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: None)
        with pytest.raises(RuntimeError, match="samtools not found"):
            rmm._run_align_and_sort(["bwa"], str(tmp_path / "o.bam"), 1, "bwa")

    def test_map_dispatch(self, monkeypatch, tmp_path):
        import gside.engine.read_mapper as rmm

        seen = {}

        class _Mapper:
            def map(self, *a, **kw):
                seen.update(kw)
                return {"ok": True}

        monkeypatch.setattr(rmm, "_get_mapper", lambda name: _Mapper())
        out = ReadMapper.map(["a.fq"], "ref", "o.bam", mode="bwa", extra=1)
        assert out == {"ok": True} and seen["extra"] == 1


class TestMapperMap:
    def _stub_env(self, monkeypatch):
        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: f"/bin/{tool}")
        monkeypatch.setattr(rmm, "_ensure_bwa_index", lambda ref: None)
        return rmm

    def test_bwa_cmd_and_result(self, monkeypatch, tmp_path):
        rmm = self._stub_env(monkeypatch)
        seen = {}
        monkeypatch.setattr(
            rmm, "_run_align_and_sort", lambda cmd, out, t, name: seen.update(cmd=cmd)
        )
        out = rmm.BwaReadMapper().map(["a.fq", "b.fq"], "ref.fna", "o.bam", extra_args="--x 1")
        assert seen["cmd"][:4] == ["/bin/bwa", "mem", "-t", seen["cmd"][3]]
        assert "--x" in seen["cmd"] and out["paired_end"] is True
        assert out["indexed"] is False

    def test_bwa_missing_binary(self, monkeypatch):
        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: None)
        with pytest.raises(RuntimeError, match="bwa not found"):
            rmm.BwaReadMapper().map(["a.fq"], "ref", "o.bam")

    def test_minimap2_cmd(self, monkeypatch):
        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: "/bin/minimap2")
        seen = {}
        monkeypatch.setattr(
            rmm, "_run_align_and_sort", lambda cmd, out, t, name: seen.update(cmd=cmd)
        )
        out = rmm.Minimap2ReadMapper().map(["a.fq"], "ref", "o.bam", preset="asm5")
        assert seen["cmd"][:4] == ["/bin/minimap2", "-ax", "asm5", "--secondary=no"]
        assert out["aligner"] == "minimap2" and out["preset"] == "asm5"

    def test_minimap2_missing_binary(self, monkeypatch):
        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: None)
        with pytest.raises(RuntimeError, match="minimap2 not found"):
            rmm.Minimap2ReadMapper().map(["a.fq"], "ref", "o.bam")


class TestRunAlignAndSort:
    def _popen(self, monkeypatch, rc=0, err=""):
        import subprocess

        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: "/bin/samtools")

        class _Ap:
            stdout = object()

            def __init__(self, *a, **k):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def wait(self, timeout=None):
                return rc

        monkeypatch.setattr(subprocess, "Popen", _Ap)
        return rmm

    def test_success(self, monkeypatch, tmp_path):
        import subprocess
        import types

        rmm = self._popen(monkeypatch)
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        rmm._run_align_and_sort(["bwa"], str(tmp_path / "o.bam"), 1, "bwa mem")

    def test_aligner_failure(self, monkeypatch, tmp_path):
        import subprocess
        import types

        rmm = self._popen(monkeypatch, rc=2)
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        with pytest.raises(RuntimeError, match="bwa mem failed"):
            rmm._run_align_and_sort(["bwa"], str(tmp_path / "o.bam"), 1, "bwa mem")

    def test_sort_failure(self, monkeypatch, tmp_path):
        import subprocess
        import types

        rmm = self._popen(monkeypatch)
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=1, stdout="", stderr="sort-broke"),
        )
        with pytest.raises(RuntimeError, match="samtools sort failed"):
            rmm._run_align_and_sort(["bwa"], str(tmp_path / "o.bam"), 1, "bwa mem")


class TestEnsureIndexBuild:
    def test_builds_when_missing(self, monkeypatch, tmp_path):
        import subprocess

        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: "/bin/bwa")
        seen = {}
        monkeypatch.setattr(subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd))
        rmm._ensure_bwa_index(str(tmp_path / "ref.fna"))
        assert seen["cmd"][:3] == ["/bin/bwa", "index", str(tmp_path / "ref.fna")]


class TestSniffBreak:
    def test_stops_after_100_records(self, tmp_path):
        rec = "@r\n" + "A" * 10 + "\n+\n" + "I" * 10 + "\n"
        p = tmp_path / "many.fq"
        p.write_text(rec * 101)
        assert _sniff_read_type([str(p)]) == "short"

    def test_auto_select(self, monkeypatch):
        import gside.engine.read_mapper as rmm

        seen = {}

        class _Mapper:
            def map(self, *a, **kw):
                return seen.update(mode="m") or {"ok": True}

        monkeypatch.setattr(rmm, "_get_mapper", lambda name: _Mapper())
        assert ReadMapper.map(["x.fna"], "ref", "o.bam") == {"ok": True}

    def test_minimap2_extra(self, monkeypatch):
        import gside.engine.read_mapper as rmm

        monkeypatch.setattr(rmm, "which", lambda tool: f"/bin/{tool}")

        class _Ap:
            stdout = object()

            def __init__(self, *a, **k):
                cmds.append(a[0])

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def wait(self, timeout=None):
                return 0

        cmds = []
        monkeypatch.setattr(subprocess, "Popen", _Ap)
        monkeypatch.setattr(subprocess, "run", lambda cmd, **kw: cmds.append(cmd) or _completed(""))
        rmm.Minimap2ReadMapper().map(["a.fq"], "ref", "o.bam", extra_args="--x 1")
        assert any("--x" in cmd for cmd in cmds)
