"""Pure unit tests for multigene_identifier internals (blast stubbed)."""

from __future__ import annotations

import subprocess
import types

import gside.analysis.multigene_identifier as mg
from gside.analysis.multigene_identifier import (
    _blast_contigs,
    _db_version,
    _load_rules,
)


def _completed(stdout: str = "", returncode: int = 0, stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


_HIT_LINE = "\t".join(["mk~~~g1~~~x", "99.0", "400", "450"])


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
    def test_short_lines_skipped(self, monkeypatch):
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: _completed(f"short\n{_HIT_LINE}\n"))
        hits = _blast_contigs("x.fna")
        assert [h["gene"] for h in hits] == ["g1"]

    def test_blastn_fallback_literal(self, monkeypatch):
        calls = []

        def fake_run(*args, **kwargs):
            calls.append(args[0])
            if args[0][0] == "sh":
                return _completed("", 1)
            return _completed(_HIT_LINE)

        monkeypatch.setattr(subprocess, "run", fake_run)
        hits = _blast_contigs("x.fna")
        assert calls[0][0] == "sh"
        assert hits and hits[0]["gene"] == "g1"


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
