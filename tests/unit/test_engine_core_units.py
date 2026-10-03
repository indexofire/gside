"""Pure unit tests for engine facade, lazy registry, remaining backend paths."""

from __future__ import annotations

import subprocess
import types
from pathlib import Path

import pytest

from gside.engine import SequenceMatcher
from gside.engine.backends import minimap2 as mm2_mod


def _completed(stdout: str = "", returncode: int = 0, stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


class TestSelectBackend:
    def test_small_defaults_blastn(self):
        assert SequenceMatcher._select_backend("x.fna", "auto") == "blastn"

    def test_prot_goes_blastp(self):
        assert SequenceMatcher._select_backend("x.faa", "prot") == "blastp"

    def test_big_file_goes_minimap2(self, tmp_path):
        big = tmp_path / "big.fna"
        with open(big, "wb") as fh:
            fh.truncate(11_000_000)
        assert SequenceMatcher._select_backend(big, "auto") == "minimap2"

    def test_missing_file_defaults_blastn(self, tmp_path):
        assert SequenceMatcher._select_backend(tmp_path / "nope.fna", "auto") == "blastn"


class TestMatchDispatch:
    def test_kmer_modes_rejected(self):
        with pytest.raises(ValueError, match="KmerDistance"):
            SequenceMatcher.match("q.fna", "db", mode="mash")

    def test_blast_dispatch(self, monkeypatch):
        import gside.engine as eng

        seen = {}

        class _Backend:
            def find(self, **kw):
                seen.update(kw)
                return ["hit"]

        monkeypatch.setattr(eng, "get_backend", lambda name, **kw: _Backend())
        assert SequenceMatcher.match("q.fna", "db", mode="blastn") == ["hit"]
        assert seen["db_path"] == "db"

    def test_minimap2_dispatch(self, monkeypatch):
        import gside.engine as eng

        seen = {}

        class _Backend:
            def find(self, **kw):
                seen.update(kw)
                return []

        monkeypatch.setattr(eng, "get_backend", lambda name, **kw: _Backend())
        SequenceMatcher.match("q.fna", "db", mode="minimap2")
        assert seen["target"] == Path("db")

    def test_blastp_tool_kwarg(self, monkeypatch):
        import gside.engine as eng

        seen = {}

        def fake_get(name, **kw):
            seen.update(kw)

            class _Backend:
                def find(self, **kw):
                    return []

            return _Backend()

        monkeypatch.setattr(eng, "get_backend", fake_get)
        SequenceMatcher.match("q.faa", "db", mode="blastp")
        assert seen.get("tool") == "blastp"


class TestLazyEnsure:
    def test_builtin_lazy_import(self, monkeypatch):
        import gside.engine.backends as backends

        monkeypatch.setattr(backends, "_BUILTINS", {"fake-x": ("pathlib", "Path")})
        assert backends.get_backend("fake-x") == Path()


class TestMinimapExtras:
    def test_make_index(self, monkeypatch, tmp_path):
        from gside.engine.backends.minimap2 import MinimapBackend

        monkeypatch.setattr(mm2_mod, "which", lambda tool: "/bin/minimap2")
        seen = {}
        monkeypatch.setattr(subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd))
        MinimapBackend().make_index(tmp_path / "r.fna", tmp_path / "r.mmi")
        assert seen["cmd"][:3] == ["/bin/minimap2", "-d", str(tmp_path / "r.mmi")]

    def test_bool_kwarg(self, monkeypatch):
        from gside.engine.backends.minimap2 import MinimapBackend

        monkeypatch.setattr(mm2_mod, "which", lambda tool: "/bin/minimap2")
        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        MinimapBackend().find(Path("q.fna"), Path("r.fna"), split_prefix=True)
        assert "-split_prefix" in seen["cmd"]

    def test_blank_and_garbage_lines(self, monkeypatch):
        from gside.engine.backends.minimap2 import MinimapBackend

        monkeypatch.setattr(mm2_mod, "which", lambda tool: "/bin/minimap2")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: _completed("\nnot-a-paf-line\n"),
        )
        assert MinimapBackend().find(Path("q.fna"), Path("r.fna")) == []


class TestMmseqsErrors:
    def test_cluster_failure(self, monkeypatch, tmp_path):
        import gside.engine.backends.mmseqs2 as mm2
        from gside.engine.backends.mmseqs2 import Mmseqs2Backend

        monkeypatch.setattr(mm2, "which", lambda tool: "/bin/mmseqs")
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("", 1, "boom"))
        with pytest.raises(RuntimeError, match="easy-linclust failed"):
            Mmseqs2Backend().cluster(tmp_path / "a.fna", tmp_path / "out")

    def test_cluster_missing_tsv(self, monkeypatch, tmp_path):
        import gside.engine.backends.mmseqs2 as mm2
        from gside.engine.backends.mmseqs2 import Mmseqs2Backend

        monkeypatch.setattr(mm2, "which", lambda tool: "/bin/mmseqs")
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(""))
        with pytest.raises(RuntimeError, match="did not produce"):
            Mmseqs2Backend().cluster(tmp_path / "a.fna", tmp_path / "out")


class TestMatchAuto:
    def test_auto_selects_and_dispatches(self, monkeypatch):
        import gside.engine as eng

        seen = {}

        class _Backend:
            def find(self, **kw):
                seen.update(kw)
                return []

        monkeypatch.setattr(eng, "get_backend", lambda name, **kw: _Backend())
        SequenceMatcher.match("q.fna", "db")
        assert seen["db_path"] == "db"
