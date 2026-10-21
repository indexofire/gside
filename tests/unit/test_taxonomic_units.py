"""Unit tests for the CheckM2 + GTDB-Tk assembly validator (binaries stubbed)."""

from __future__ import annotations

import json
import types
from datetime import datetime
from pathlib import Path

import pytest

import gside.analysis.taxonomic_validator as tv
from gside import __version__
from gside.analysis.taxonomic_validator import (
    SCHEMA,
    ValidateResult,
    _build_interpretation,
    _quality_verdict,
    _run_checkm2,
    _run_gtdbtk,
    _tool_version,
    validate_genome,
)


def _completed(returncode: int = 0, stdout: str = "", stderr: str = ""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


class TestToolVersion:
    def test_first_line(self, monkeypatch):
        monkeypatch.setattr(
            tv.subprocess, "run", lambda *a, **k: _completed(stdout="1.0.2\nextra\n")
        )
        assert _tool_version("/bin/checkm2") == "1.0.2"

    def test_falls_back_to_stderr(self, monkeypatch):
        monkeypatch.setattr(tv.subprocess, "run", lambda *a, **k: _completed(stderr="v2.4.0"))
        assert _tool_version("/bin/gtdbtk") == "v2.4.0"

    def test_missing_binary_is_empty(self, monkeypatch):
        def boom(*args, **kwargs):
            raise OSError("missing")

        monkeypatch.setattr(tv.subprocess, "run", boom)
        assert _tool_version("/bin/nope") == ""


class TestCheckm2:
    def _ready(self, monkeypatch, tmp_path):
        monkeypatch.setattr(tv, "CHECKM2_DB", tmp_path / "db")
        monkeypatch.setattr(tv, "_find_tool", lambda name: "/bin/checkm2")
        monkeypatch.setattr(tv, "_tool_version", lambda binary: "1.0.2")

    def test_missing_database(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "CHECKM2_DB", None)
        result = _run_checkm2("x.fna", tmp_path)
        assert result["status"] == "missing"
        assert "CHECKM2DB" in result["error"]

    def test_missing_binary(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "CHECKM2_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: None)
        result = _run_checkm2("x.fna", tmp_path)
        assert result["status"] == "missing"
        assert "not found" in result["error"]

    def test_success(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)

        def fake(cmd, **kwargs):
            out = Path(cmd[cmd.index("--output") + 1])
            out.write_text("Name\tCompleteness\tContamination\ng1\t99.0\t1.0\n")
            return _completed()

        monkeypatch.setattr(tv.subprocess, "run", fake)
        result = _run_checkm2("x.fna", tmp_path / "o", threads=8)
        assert result["status"] == "ok"
        assert result["completeness"] == 99.0
        assert result["contamination"] == 1.0
        assert result["version"] == "1.0.2"

    def test_nonzero_exit(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)
        monkeypatch.setattr(
            tv.subprocess, "run", lambda *a, **k: _completed(returncode=2, stderr="boom")
        )
        result = _run_checkm2("x.fna", tmp_path)
        assert result["status"] == "failed"
        assert "exit 2" in result["error"]

    def test_missing_report(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)
        monkeypatch.setattr(tv.subprocess, "run", lambda *a, **k: _completed())
        result = _run_checkm2("x.fna", tmp_path)
        assert result["status"] == "failed"
        assert "no checkm2_results.tsv" in result["error"]

    def test_non_numeric_report(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)

        def fake(cmd, **kwargs):
            out = Path(cmd[cmd.index("--output") + 1])
            out.write_text("Name\tCompleteness\tContamination\ng1\tbad\tbad\n")
            return _completed()

        monkeypatch.setattr(tv.subprocess, "run", fake)
        result = _run_checkm2("x.fna", tmp_path)
        assert result["status"] == "failed"
        assert "numeric" in result["error"]


class TestGtdbtk:
    def _ready(self, monkeypatch, tmp_path):
        monkeypatch.setattr(tv, "GTDB_DB", tmp_path / "gtdb")
        monkeypatch.setattr(tv, "_find_tool", lambda name: "/bin/gtdbtk")
        monkeypatch.setattr(tv, "_tool_version", lambda binary: "2.4.0")

    def test_missing_database(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "GTDB_DB", None)
        result = _run_gtdbtk("x.fna", tmp_path)
        assert result["status"] == "missing"
        assert "GTDBTK_DATA_PATH" in result["note"]

    def test_missing_binary(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "GTDB_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: None)
        result = _run_gtdbtk("x.fna", tmp_path)
        assert result["status"] == "missing"
        assert "not found" in result["note"]

    def test_bac120_success(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)
        (tmp_path / "q.fna").write_text(">q\nACGT\n")

        def fake(cmd, **kwargs):
            out = Path(cmd[cmd.index("--out_dir") + 1])
            out.mkdir(parents=True, exist_ok=True)
            (out / "gtdbtk.bac120.summary.tsv").write_text(
                "user_genome\tclassification\nq\td__Bacteria;s__X\n"
            )
            return _completed()

        monkeypatch.setattr(tv.subprocess, "run", fake)
        result = _run_gtdbtk(str(tmp_path / "q.fna"), tmp_path / "o", threads=8)
        assert result["status"] == "ok"
        assert result["classification"] == "d__Bacteria;s__X"

    def test_ar122_fallback(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)
        (tmp_path / "q.fna").write_text(">q\nACGT\n")

        def fake(cmd, **kwargs):
            out = Path(cmd[cmd.index("--out_dir") + 1])
            out.mkdir(parents=True, exist_ok=True)
            (out / "gtdbtk.ar122.summary.tsv").write_text(
                "user_genome\tclassification\nq\td__Archaea;s__Y\n"
            )
            return _completed()

        monkeypatch.setattr(tv.subprocess, "run", fake)
        assert _run_gtdbtk(str(tmp_path / "q.fna"), tmp_path / "o")["classification"] == (
            "d__Archaea;s__Y"
        )

    def test_nonzero_exit(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)
        (tmp_path / "q.fna").write_text(">q\nACGT\n")
        monkeypatch.setattr(
            tv.subprocess, "run", lambda *a, **k: _completed(returncode=1, stderr="boom")
        )
        result = _run_gtdbtk(str(tmp_path / "q.fna"), tmp_path / "o")
        assert result["status"] == "failed"
        assert "exit 1" in result["note"]

    def test_missing_summary(self, tmp_path, monkeypatch):
        self._ready(monkeypatch, tmp_path)
        (tmp_path / "q.fna").write_text(">q\nACGT\n")
        monkeypatch.setattr(tv.subprocess, "run", lambda *a, **k: _completed())
        result = _run_gtdbtk(str(tmp_path / "q.fna"), tmp_path / "o")
        assert result["status"] == "failed"
        assert "no bac120/ar122 summary" in result["note"]


class TestQualityVerdict:
    @pytest.mark.parametrize(
        "completeness,contamination,verdict",
        [
            (99.0, 1.0, "PASS"),
            (95.0, 5.0, "PASS"),
            (30.0, 1.0, "WARNING"),
            (99.0, 12.0, "WARNING"),
            (80.0, 1.0, "CHECK"),
            (None, 1.0, "UNKNOWN"),
        ],
    )
    def test_bands(self, completeness, contamination, verdict):
        assert _quality_verdict(completeness, contamination) == verdict


class TestInterpretation:
    def test_pass_with_taxonomy(self):
        result = ValidateResult(
            quality={"verdict": "PASS", "completeness": 99.0, "contamination": 1.0},
            taxonomy={"classification": "d__Bacteria;s__X", "note": ""},
        )
        text = _build_interpretation(result)
        assert "Genome quality: PASS" in text
        assert "GTDB-Tk: d__Bacteria;s__X" in text

    def test_missing_checkm2(self):
        result = ValidateResult(
            quality={"verdict": "UNKNOWN", "completeness": None, "contamination": None},
            taxonomy={"classification": "", "note": ""},
            status="partial",
        )
        text = _build_interpretation(result)
        assert "CheckM2 unavailable" in text
        assert "GTDB-Tk: unavailable" in text
        assert "Status: partial" in text


class TestValidateGenome:
    def _stub(self, monkeypatch, checkm2, gtdbtk):
        monkeypatch.setattr(tv, "CHECKM2_DB", Path("/db/checkm2"))
        monkeypatch.setattr(tv, "GTDB_DB", Path("/db/gtdb"))
        monkeypatch.setattr(tv, "_run_checkm2", lambda c, o, t: checkm2)
        monkeypatch.setattr(tv, "_run_gtdbtk", lambda c, o, t: gtdbtk)

    def test_ok_writes_full_report(self, tmp_path, monkeypatch):
        self._stub(
            monkeypatch,
            {
                "status": "ok",
                "completeness": 99.0,
                "contamination": 1.0,
                "version": "1.0.2",
                "error": "",
            },
            {"status": "ok", "classification": "d__Bacteria;s__X", "version": "2.4.0", "note": ""},
        )
        source = tmp_path / "GCF_1_genomic.fna"
        source.write_text(">c\nACGT\n")
        out = tmp_path / "tax"
        result = validate_genome(source, out, threads=8)

        assert result.status == "ok"
        assert result.genome == {"name": "GCF_1_genomic", "path": str(source.resolve())}
        assert result.parameters["threads"] == 8
        assert result.parameters["checkm2_db"] == "/db/checkm2"
        assert result.software == {"gside": __version__, "checkm2": "1.0.2", "gtdbtk": "2.4.0"}
        assert result.tools == {"checkm2": "ok", "gtdbtk": "ok"}
        assert result.quality == {"completeness": 99.0, "contamination": 1.0, "verdict": "PASS"}
        assert result.taxonomy["classification"] == "d__Bacteria;s__X"
        assert result.errors == []
        datetime.fromisoformat(result.created_at)

        payload = json.loads((out / "validation.json").read_text())
        assert payload["schema"] == SCHEMA
        assert payload["genome"]["name"] == "GCF_1_genomic"
        assert payload["quality"]["verdict"] == "PASS"

    def test_partial_when_gtdbtk_missing(self, tmp_path, monkeypatch):
        self._stub(
            monkeypatch,
            {
                "status": "ok",
                "completeness": 99.0,
                "contamination": 1.0,
                "version": "1.0.2",
                "error": "",
            },
            {"status": "missing", "classification": "", "version": "", "note": "no gtdb"},
        )
        source = tmp_path / "g.fna"
        source.write_text(">c\nACGT\n")
        result = validate_genome(source, tmp_path / "tax")
        assert result.status == "partial"
        assert result.tools["gtdbtk"] == "missing"
        assert result.errors == ["no gtdb"]
        assert result.taxonomy["note"] == "no gtdb"

    def test_error_when_both_missing(self, tmp_path, monkeypatch):
        self._stub(
            monkeypatch,
            {
                "status": "missing",
                "completeness": None,
                "contamination": None,
                "version": "",
                "error": "no checkm2db",
            },
            {"status": "missing", "classification": "", "version": "", "note": "no gtdb"},
        )
        source = tmp_path / "g.fna"
        source.write_text(">c\nACGT\n")
        result = validate_genome(source, tmp_path / "nested" / "tax")
        assert result.status == "error"
        assert result.quality["verdict"] == "UNKNOWN"
        assert len(result.errors) == 2
        assert (tmp_path / "nested" / "tax" / "validation.json").is_file()


class TestFindTool:
    def test_resolves_via_path(self, monkeypatch):
        monkeypatch.setattr(
            tv.shutil,
            "which",
            lambda name, path=None: f"/bin/{name}" if name == "checkm2" else None,
        )
        assert tv._find_tool("checkm2") == "/bin/checkm2"
        assert tv._find_tool("gtdbtk") is None


class TestRunnerFailures:
    def test_checkm2_spawn_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "CHECKM2_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: "/bin/checkm2")
        monkeypatch.setattr(tv, "_tool_version", lambda binary: "")

        def boom(*args, **kwargs):
            raise OSError("cannot exec")

        monkeypatch.setattr(tv.subprocess, "run", boom)
        result = tv._run_checkm2("x.fna", tmp_path)
        assert result["status"] == "failed"
        assert "cannot exec" in result["error"]

    def test_checkm2_header_only_report(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "CHECKM2_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: "/bin/checkm2")
        monkeypatch.setattr(tv, "_tool_version", lambda binary: "")

        def fake(cmd, **kwargs):
            Path(cmd[cmd.index("--output") + 1]).write_text("Name\tCompleteness\n")
            return _completed()

        monkeypatch.setattr(tv.subprocess, "run", fake)
        result = tv._run_checkm2("x.fna", tmp_path)
        assert result["status"] == "failed"
        assert "numeric" in result["error"]

    def test_gtdbtk_spawn_error(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "GTDB_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: "/bin/gtdbtk")
        monkeypatch.setattr(tv, "_tool_version", lambda binary: "")
        (tmp_path / "q.fna").write_text(">q\nACGT\n")

        def boom(*args, **kwargs):
            raise OSError("cannot exec")

        monkeypatch.setattr(tv.subprocess, "run", boom)
        result = tv._run_gtdbtk(str(tmp_path / "q.fna"), tmp_path / "o")
        assert result["status"] == "failed"
        assert "cannot exec" in result["note"]

    def test_gtdbtk_blank_classification(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "GTDB_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: "/bin/gtdbtk")
        monkeypatch.setattr(tv, "_tool_version", lambda binary: "")
        (tmp_path / "q.fna").write_text(">q\nACGT\n")

        def fake(cmd, **kwargs):
            out = Path(cmd[cmd.index("--out_dir") + 1])
            out.mkdir(parents=True, exist_ok=True)
            (out / "gtdbtk.bac120.summary.tsv").write_text("user_genome\tclassification\nq\t\n")
            return _completed()

        monkeypatch.setattr(tv.subprocess, "run", fake)
        result = tv._run_gtdbtk(str(tmp_path / "q.fna"), tmp_path / "o")
        assert result["status"] == "failed"
        assert "no classification" in result["note"]
