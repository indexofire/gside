"""Pure unit tests for skani/blast backends (subprocess stubbed)."""

from __future__ import annotations

import subprocess
import types
from pathlib import Path

import pytest

from gside.engine.backends import blast as blast_mod
from gside.engine.backends import skani as skani_mod
from gside.engine.backends.blast import BlastBackend
from gside.engine.backends.skani import SkaniBackend


def _completed(stdout: str = "", returncode: int = 0, stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


@pytest.fixture()
def nobin(monkeypatch):
    monkeypatch.setattr(skani_mod, "which", lambda tool: None)
    monkeypatch.setattr(blast_mod, "which", lambda tool: None)


@pytest.fixture()
def fakebin(monkeypatch):
    monkeypatch.setattr(skani_mod, "which", lambda tool: "/bin/skani")
    monkeypatch.setattr(blast_mod, "which", lambda tool: "/bin/blastn")
    import gside.engine.backends.kmer as kmer_mod

    monkeypatch.setattr(kmer_mod, "which", lambda tool: "/bin/mash")


class TestFindBin:
    def test_skani_missing(self, nobin):
        with pytest.raises(RuntimeError, match="skani not found"):
            SkaniBackend()

    def test_blast_missing(self, nobin):
        with pytest.raises(RuntimeError, match="blastn not found"):
            BlastBackend()


class TestExplicitBinary:
    def test_blast_binary_bypasses_discovery(self, nobin, monkeypatch):
        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        backend = BlastBackend(binary="/explicit/blastn")
        backend.find(Path("q.fna"), "db")
        assert seen["cmd"][0] == "/explicit/blastn"


def _skani_table(header: str, rows: list[str]) -> str:
    return "\n".join([header] + rows)


_NAMED_HEADER = "Ref_file\tQuery_file\tANI\tAlign_fraction_query\tfoo"
_LEGACY_HEADER = "a\tb\tc\td\te\tf"


class TestSkaniSearch:
    def test_named_header(self, fakebin, monkeypatch):
        table = _skani_table(_NAMED_HEADER, ["g1\tq\t99.0\t0.9"])
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(table))
        hits = SkaniBackend().search(Path("q.fna"), Path("db"))
        assert [(h.ref, h.ani, h.aligned_fraction) for h in hits] == [("g1", 99.0, 0.9)]

    def test_legacy_fallback_indices(self, fakebin, monkeypatch):
        table = _skani_table(_LEGACY_HEADER, ["g1\tx\ty\t0.8\tz\t97.5"])
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(table))
        hits = SkaniBackend().search(Path("q.fna"), Path("db"))
        assert [(h.ref, h.ani, h.aligned_fraction) for h in hits] == [("g1", 97.5, 0.8)]

    def test_empty_and_short_rows(self, fakebin, monkeypatch):
        table = _skani_table(_NAMED_HEADER, ["", "only\two", "g2\tq\t98.0\t0.7"])
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(table))
        hits = SkaniBackend().search(Path("q.fna"), Path("db"))
        assert [h.ref for h in hits] == ["g2"]

    def test_no_output(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("\n"))
        assert SkaniBackend().search(Path("q.fna"), Path("db")) == []

    def test_failure_raises(self, fakebin, monkeypatch):
        monkeypatch.setattr(
            subprocess, "run", lambda *a, **k: _completed("", returncode=1, stderr="boom")
        )
        with pytest.raises(RuntimeError, match="skani search failed"):
            SkaniBackend().search(Path("q.fna"), Path("db"))


_BLAST_LINE = "\t".join(
    ["q1", "s1", "99.0", "400", "2", "0", "10", "409", "1", "400", "1e-50", "700", "500", "450"]
)


class TestBlastFind:
    def _backend(self, fakebin, monkeypatch, stdout=_BLAST_LINE, returncode=0):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(stdout, returncode))
        return BlastBackend()

    def test_parse_and_filter(self, fakebin, monkeypatch):
        hits = self._backend(fakebin, monkeypatch).find(
            Path("q.fna"), "db", min_identity=90.0, min_coverage=50.0
        )
        assert len(hits) == 1 and hits[0].query_id == "q1"

    def test_filter_rejects(self, fakebin, monkeypatch):
        hits = self._backend(fakebin, monkeypatch).find(Path("q.fna"), "db", min_identity=99.9)
        assert hits == []

    def test_kwargs_mapping(self, fakebin, monkeypatch):
        seen = {}

        def fake_run(cmd, **kw):
            seen["cmd"] = cmd
            return _completed("")

        monkeypatch.setattr(subprocess, "run", fake_run)
        BlastBackend().find(
            Path("q.fna"),
            "db",
            pident=97,
            num_threads=8,
            threads=2,
            dust="no",
            keep=None,
        )
        cmd = seen["cmd"]
        assert "-perc_identity" in cmd and "97" in cmd
        assert cmd[cmd.index("-num_threads") + 1] == "8"
        assert "-dust" in cmd and "no" in cmd
        assert "keep" not in " ".join(cmd)

    def test_failure_raises(self, fakebin, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("", 1, "boom"))
        with pytest.raises(RuntimeError, match="exit 1"):
            BlastBackend().find(Path("q.fna"), "db")


class TestEnsureIndex:
    def test_existing_index_noop(self, fakebin, tmp_path, monkeypatch):
        (tmp_path / "db.nhr").write_bytes(b"x")
        called = []
        monkeypatch.setattr(BlastBackend, "make_db", lambda self, f, p, t="nucl": called.append(f))
        BlastBackend().ensure_index(str(tmp_path / "db"))
        assert called == []

    def test_builds_from_candidate(self, fakebin, tmp_path, monkeypatch):
        (tmp_path / "db.fna").write_text(">a\nACGT\n")
        called = []
        monkeypatch.setattr(
            BlastBackend, "make_db", lambda self, f, p, t="nucl": called.append((f, p))
        )
        BlastBackend().ensure_index(str(tmp_path / "db"))
        assert called[0][0] == tmp_path / "db.fna"

    def test_no_source_raises(self, fakebin, tmp_path):
        with pytest.raises(FileNotFoundError):
            BlastBackend().ensure_index(str(tmp_path / "db"))

    def test_make_db_missing_binary(self, nobin):
        backend = BlastBackend.__new__(BlastBackend)
        with pytest.raises(RuntimeError, match="makeblastdb not found"):
            backend.make_db(Path("a"), Path("b"))


class TestMakeDbAndFlags:
    def test_make_db_cmd(self, fakebin, monkeypatch, tmp_path):
        import subprocess

        seen = {}
        monkeypatch.setattr(subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd))
        BlastBackend().make_db(tmp_path / "a.fna", tmp_path / "db")
        assert seen["cmd"][:4] == ["/bin/blastn", "-in", str(tmp_path / "a.fna"), "-dbtype"]

    def test_bool_flag(self, fakebin, monkeypatch):
        import subprocess

        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        BlastBackend().find(Path("q.fna"), "db", dust=True)
        assert "-dust" in seen["cmd"]

    def test_blank_and_garbage_lines(self, fakebin, monkeypatch):
        import subprocess

        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: _completed("\nnot-enough-cols\n" + _BLAST_LINE),
        )
        hits = BlastBackend().find(Path("q.fna"), "db")
        assert len(hits) == 1


class TestBackendThreads:
    def test_skani_search_carries_t_flag(self, fakebin, monkeypatch):
        from gside.engine.backends.skani import SkaniBackend

        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("\n")
        )
        SkaniBackend(threads=3).search(Path("q.fna"), Path("db"))
        i = seen["cmd"].index("-t")
        assert seen["cmd"][i + 1] == "3"

    def test_mash_distance_default_threads_flag(self, fakebin, monkeypatch):
        from gside.engine.backends.kmer import MashBackend

        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        MashBackend().distance(Path("q.msh"), Path("r.msh"))
        i = seen["cmd"].index("-p")
        assert seen["cmd"][i + 1] == "4"

    def test_mash_explicit_p_kwarg_overrides_threads(self, fakebin, monkeypatch):
        from gside.engine.backends.kmer import MashBackend

        seen = {}
        monkeypatch.setattr(
            subprocess, "run", lambda cmd, **kw: seen.update(cmd=cmd) or _completed("")
        )
        MashBackend(threads=6).distance(Path("q.msh"), Path("r.msh"), p=500)
        flags = [(f, seen["cmd"][i + 1]) for i, f in enumerate(seen["cmd"]) if f == "-p"]
        assert flags == [("-p", "500")]
