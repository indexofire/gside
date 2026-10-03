"""Pure unit tests for kmer/minimap2 backends (subprocess stubbed)."""

from __future__ import annotations

import subprocess
import types
from pathlib import Path

import pytest

from gside.engine.backends import kmer as kmer_mod
from gside.engine.backends import minimap2 as mm2_mod
from gside.engine.backends.kmer import MashBackend, SourmashBackend
from gside.engine.backends.minimap2 import MinimapBackend


def _completed(stdout: str = "", returncode: int = 0, stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


@pytest.fixture()
def fakebin(monkeypatch):
    monkeypatch.setattr(kmer_mod, "which", lambda tool: f"/bin/{tool}")
    monkeypatch.setattr(mm2_mod, "which", lambda tool: "/bin/minimap2")


@pytest.fixture()
def nobin(monkeypatch):
    monkeypatch.setattr(kmer_mod, "which", lambda tool: None)
    monkeypatch.setattr(mm2_mod, "which", lambda tool: None)


class TestFindBin:
    def test_mash_missing(self, nobin):
        with pytest.raises(RuntimeError, match="mash not found"):
            MashBackend()

    def test_sourmash_missing(self, nobin):
        with pytest.raises(RuntimeError, match="sourmash not found"):
            SourmashBackend()

    def test_minimap2_missing(self, nobin):
        with pytest.raises(RuntimeError, match="minimap2 not found"):
            MinimapBackend()


_MASH_LINE = "refA\tquery1\t0.02\t0.0\t80/1000"


class TestMashDistance:
    def test_parse_slash_hashes(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(_MASH_LINE))
        (res,) = MashBackend().distance(Path("q.msh"), Path("r.msh"))
        assert (res.reference_id, res.query_id) == ("refA", "query1")
        assert res.distance == 0.02 and res.shared_hashes == 80 and res.total_hashes == 1000

    def test_parse_split_hashes(self, fakebin, monkeypatch):
        line = "refB\tquery1\t0.05\t0.001\t40\t1000\textra"
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(line))
        (res,) = MashBackend().distance(Path("q.msh"), Path("r.msh"))
        assert (res.shared_hashes, res.total_hashes) == (40, 1000)

    def test_skips_comments_shorts_garbage(self, fakebin, monkeypatch):
        stdout = "# comment\n\nshort\tline\nrefC\tq\tbadfloat\t0.0\t1/10\n" + _MASH_LINE
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(stdout))
        res = MashBackend().distance(Path("q.msh"), Path("r.msh"))
        assert [r.reference_id for r in res] == ["refA"]

    def test_failure_raises(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("", 1, "boom"))
        with pytest.raises(RuntimeError, match="mash dist failed"):
            MashBackend().distance(Path("q.msh"), Path("r.msh"))

    def test_kwargs_appended(self, fakebin, monkeypatch):
        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        MashBackend().distance(Path("q.msh"), Path("r.msh"), p=500)
        assert "-p" in seen["cmd"] and "500" in seen["cmd"]

    def test_sketch_cmd(self, fakebin, monkeypatch, tmp_path):
        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        out = MashBackend().sketch(tmp_path / "g.fna", tmp_path / "g", individual=True)
        assert seen["cmd"][:6] == ["/bin/mash", "sketch", "-k", "21", "-s", "1000"]
        assert "-i" in seen["cmd"]
        assert out == tmp_path / "g.msh"


_SMASH_CSV = (
    "name,query_name,intersect_bp,containment,extra1,extra2,extra3\n"
    "q1,ok,fileX,refX,x,0.9,y\n"
    "short,row\n"
    "refY,q1,x,notanum,y\n"
)


class TestSourmashDistance:
    def test_parse(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(_SMASH_CSV))
        res = SourmashBackend().distance(Path("q.sig"), Path("r.sig"))
        assert [(r.reference_id, r.distance) for r in res] == [("refX", 0.1)]
        assert res[0].backend == "sourmash"

    def test_failure_raises(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("", 1, "boom"))
        with pytest.raises(RuntimeError, match="sourmash search failed"):
            SourmashBackend().distance(Path("q.sig"), Path("r.sig"))

    def test_sketch_name_flag(self, fakebin, monkeypatch, tmp_path):
        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        SourmashBackend().sketch(tmp_path / "g.fna", tmp_path / "g.sig", name="n1")
        assert "--name" in seen["cmd"] and "n1" in seen["cmd"]


_PAF = "\t".join(
    ["q1", "1000", "100", "600", "+", "ref1", "2000", "300", "800", "490", "500", "60"]
)


class TestMinimapFind:
    def test_parse_and_filter(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(_PAF))
        hits = MinimapBackend().find(
            Path("q.fna"), Path("r.fna"), min_identity=90.0, min_coverage=40.0
        )
        assert len(hits) == 1 and hits[0].subject_id == "ref1"

    def test_filter_rejects(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(_PAF))
        assert MinimapBackend().find(Path("q.fna"), Path("r.fna"), min_identity=99.9) == []

    def test_kwargs_mapping(self, fakebin, monkeypatch):
        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        MinimapBackend().find(Path("q.fna"), Path("r.fna"), kmer=19, threads=8, drop=None)
        cmd = seen["cmd"]
        assert "-k" in cmd and "19" in cmd
        assert cmd[cmd.index("-t") + 1] == "8"

    def test_failure_raises(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("", 1, "boom"))
        with pytest.raises(RuntimeError, match="minimap2 failed"):
            MinimapBackend().find(Path("q.fna"), Path("r.fna"))


class TestScreenDelegate:
    def test_screen_calls_distance(self, fakebin, monkeypatch):
        seen = {}

        class _M(MashBackend):
            def distance(self, query, reference, max_distance=0.1, **kw):
                seen.update(max_distance=max_distance)
                return ["d"]

        monkeypatch.setattr("gside.engine.backends.kmer.MashBackend", _M)
        assert _M().screen("q", "r") == ["d"]
        assert seen["max_distance"] == 0.1


class TestSourmashShortRows:
    def test_headerless_and_short(self, fakebin, monkeypatch):
        import subprocess

        stdout = "q1,ok,fileX,refX,x,0.9,y\nshort,row\n"
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(stdout))
        res = SourmashBackend().distance(Path("q.sig"), Path("r.sig"))
        assert [(r.reference_id, r.distance) for r in res] == [("refX", 0.1)]
