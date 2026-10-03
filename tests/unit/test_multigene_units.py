"""Pure unit tests for multigene_identifier internals (blast stubbed)."""

from __future__ import annotations

import subprocess
import types
from pathlib import Path

import pytest

import gside.analysis.multigene_identifier as mg
from gside import config
from gside.analysis.multigene_identifier import (
    _blast_contigs,
    _db_version,
    _load_rules,
)


def _completed(stdout: str = "", returncode: int = 0, stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


def _blastn_stub(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, script: str) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    stub = bin_dir / "blastn"
    stub.write_text(script)
    stub.chmod(0o755)
    monkeypatch.setattr(config, "pixi_path", lambda: str(bin_dir))
    return stub


_HIT_LINE = "\t".join(
    [
        "q1",
        "mk~~~g1~~~x",
        "99.0",
        "400",
        "2",
        "0",
        "10",
        "409",
        "1",
        "400",
        "1e-50",
        "700",
        "500",
        "450",
    ]
)


class TestDbVersion:
    def test_missing_fasta_is_unknown(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mg, "_MARKERS_FASTA", tmp_path / "nope.fasta")
        assert _db_version() == "unknown"


class TestLoadRules:
    def test_missing_file_is_empty(self, monkeypatch, tmp_path):
        monkeypatch.setattr(mg, "_MARKER_RULES", tmp_path / "nope.yaml")
        assert _load_rules() == []


class TestPriorityDeadCode:
    def test_priority_block_runs_without_effect(self, monkeypatch):
        monkeypatch.setattr(
            mg,
            "_load_rules",
            lambda: [
                {
                    "species": "X",
                    "genes": ["g1"],
                    "min_hits": 1,
                    "priority_over": ["Y"],
                }
            ],
        )
        monkeypatch.setattr(
            mg,
            "_blast_contigs",
            lambda contigs: [{"gene": "g1", "identity": 99.0, "coverage": 100.0}],
        )
        from gside.analysis.multigene_identifier import identify_multigene

        assert identify_multigene("x.fna").to_dict()["species"] == "X"


class TestBlastParsing:
    def test_short_lines_skipped(self, monkeypatch, tmp_path):
        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n")
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(f"short\n{_HIT_LINE}\n"))
        hits = _blast_contigs("x.fna")
        assert [h["gene"] for h in hits] == ["g1"]


class TestBlastFailureSurface:
    def test_nonzero_exit_raises_runtime_error(self, tmp_path, monkeypatch):
        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\necho 'BLAST engine error' >&2\nexit 3\n")
        with pytest.raises(RuntimeError, match=r"blastn failed \(exit 3\)"):
            _blast_contigs("x.fna")

    def test_failure_message_carries_stderr(self, tmp_path, monkeypatch):
        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\necho 'query file corrupt' >&2\nexit 2\n")
        with pytest.raises(RuntimeError, match="query file corrupt"):
            _blast_contigs("x.fna")

    def test_discovery_failure_raises_with_fix_hint(self, monkeypatch, tmp_path):
        monkeypatch.setattr(config, "pixi_path", lambda: str(tmp_path))
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed("", 0))
        with pytest.raises(RuntimeError, match=r"blastn not found.*GSIDE_PIXI_BIN"):
            _blast_contigs("x.fna")

    def test_stub_hits_still_parsed(self, tmp_path, monkeypatch):
        _blastn_stub(
            tmp_path,
            monkeypatch,
            "#!/bin/sh\nprintf 'q1\\tmk~~~g1~~~x\\t99.0\\t400\\t2\\t0\\t10\\t409"
            "\\t1\\t400\\t1e-50\\t700\\t500\\t450\\n'\nexit 0\n",
        )
        hits = _blast_contigs("x.fna")
        assert [h["gene"] for h in hits] == ["g1"]


class TestBackendDelegation:
    def test_routes_through_blast_backend(self, tmp_path, monkeypatch):
        from gside.engine.backends import blast as blast_mod

        _blastn_stub(tmp_path, monkeypatch, "#!/bin/sh\nexit 0\n")
        calls: dict[str, object] = {}

        class _B:
            def __init__(self, tool="blastn", threads=4, binary=None):
                calls["ctor"] = binary

            def find(self, query, db_path, **kw):
                calls["find"] = (query, db_path, kw)
                return []

        monkeypatch.setattr(blast_mod, "BlastBackend", _B)
        _blast_contigs("x.fna")
        assert calls["ctor"] == str(tmp_path / "bin" / "blastn")
        query, db_path, kw = calls["find"]
        assert query == Path("x.fna")
        assert str(db_path).endswith("markers_blastdb")
        assert kw["evalue"] == 1e-10
        assert kw["word_size"] == 11


class TestMainOutput:
    def test_json_and_text(self, monkeypatch, capsys, tmp_path):
        import gside.analysis.multigene_identifier as mg

        monkeypatch.setattr(
            mg,
            "identify_multigene",
            lambda c: mg.MultiGeneResult(species="Salmonella", confidence="high"),
        )
        q = tmp_path / "q.fna"
        q.write_text(">q\nACGT\n")
        monkeypatch.setattr("sys.argv", ["prog", str(q), "--json"])
        mg.main()
        assert '"Salmonella"' in capsys.readouterr().out
        monkeypatch.setattr("sys.argv", ["prog", str(q)])
        mg.main()
        assert "Species: Salmonella (high)" in capsys.readouterr().out
