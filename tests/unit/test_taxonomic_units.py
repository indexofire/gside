"""Unit tests for taxonomic_validator (binaries stubbed, tiny real query)."""

from __future__ import annotations

import types

import gside.analysis.taxonomic_validator as tv
from gside.analysis.taxonomic_validator import (
    TaxonomyResult,
    _build_interpretation,
    _find_tool,
    _run_checkm2,
    _run_gtdbtk,
    validate_genome,
)


def _result(**kw):
    base = {
        "mode": "standard",
        "marker_gene_species": "",
        "marker_gene_confidence": "",
        "marker_gene_markers": [],
        "completeness": None,
        "contamination": None,
        "gtdb_taxonomy": "",
        "gtdb_note": "",
        "interpretation": "",
    }
    base.update(kw)
    return TaxonomyResult(**base)


class TestToDictSummary:
    def test_shape(self):
        d = _result(marker_gene_species="Salmonella").to_dict()
        assert d["marker_gene_species"] == "Salmonella"

    def test_summary_parts(self):
        s = _result(
            marker_gene_species="X",
            completeness=99.0,
            contamination=1.0,
            gtdb_taxonomy="d__Bacteria;s__X",
        ).summary
        assert "Mode: standard" in s and "Completeness: 99.0%" in s
        assert "GTDB-Tk: d__Bacteria;s__X" in s


class TestFindTool:
    def test_missing(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda *a, **k: None)
        assert _find_tool("checkm2") is None


class TestRunnersDegrade:
    def test_checkm2_no_db(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "CHECKM2_DB", None)
        assert _run_checkm2("x.fna", tmp_path) == (None, None)

    def test_checkm2_no_binary(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "CHECKM2_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: None)
        assert _run_checkm2("x.fna", tmp_path) == (None, None)

    def test_gtdbtk_no_db(self, tmp_path):
        assert _run_gtdbtk("x.fna", tmp_path) == ""

    def test_gtdbtk_no_binary(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "GTDB_DB", tmp_path)
        monkeypatch.setattr(tv, "_find_tool", lambda name: None)
        assert _run_gtdbtk("x.fna", tmp_path) == ""


class TestInterpretation:
    def test_pass(self):
        assert "PASS" in _build_interpretation(_result(completeness=99.0, contamination=1.0))

    def test_warning(self):
        assert "WARNING" in _build_interpretation(_result(completeness=30.0, contamination=1.0))

    def test_check(self):
        assert "CHECK" in _build_interpretation(_result(completeness=80.0, contamination=10.0))

    def test_agree(self):
        tax = "d__Bacteria;Salmonella enterica"
        assert "agree" in _build_interpretation(
            _result(marker_gene_species="Salmonella", gtdb_taxonomy=tax)
        )

    def test_discrepancy(self):
        assert "DISCREPANCY" in _build_interpretation(
            _result(marker_gene_species="Salmonella", gtdb_taxonomy="d__Bacteria;E. coli")
        )

    def test_empty(self):
        assert _build_interpretation(_result()) == (
            "Standard mode requested but CheckM2/GTDB-Tk databases not available"
        )


class TestValidateSimple:
    def test_simple_with_real_blastn(self, tmp_path):
        from gside.analysis.multigene_identifier import _MARKERS_FASTA

        blocks = _MARKERS_FASTA.read_text().split(">")
        inva = next(b for b in blocks if b.startswith("markers~~~inva~~~"))
        seq = inva.split("\n", 1)[1].replace("\n", "")
        q = tmp_path / "inva.fna"
        q.write_text(f">inva_q\n{seq}\n")
        out = tmp_path / "tax"
        r = validate_genome(str(q), mode="simple", output_dir=out)
        assert r.marker_gene_species == "Salmonella"
        assert "Marker gene: Salmonella" in r.interpretation
        assert (out / "validation.json").is_file()

    def test_standard_with_stubbed_runners(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tv, "_run_checkm2", lambda c, o: (99.0, 1.0))
        monkeypatch.setattr(tv, "_run_gtdbtk", lambda c, o: "d__Bacteria;Salmonella")
        import gside.analysis.multigene_identifier as mg

        monkeypatch.setattr(
            mg,
            "identify_multigene",
            lambda c: types.SimpleNamespace(
                species="Salmonella", confidence="high", detected_markers=[]
            ),
        )
        r = validate_genome(str(tmp_path / "q.fna"), mode="standard", output_dir=tmp_path / "t")
        assert "PASS" in r.interpretation and "agree" in r.interpretation


class TestValidateCommand:
    def test_help(self):
        from click.testing import CliRunner

        from gside.cli import main

        result = CliRunner().invoke(main, ["validate", "--help"])
        assert result.exit_code == 0
        assert "CheckM2" in result.output


class TestTsvParsing:
    def _checkm2_tsv(self, out, rows):
        out.mkdir(parents=True, exist_ok=True)
        (out / "checkm2_results.tsv").write_text(
            "Name\tCompleteness\tContamination\n" + "\n".join(rows) + "\n"
        )

    def test_checkm2_good(self, tmp_path, monkeypatch):
        import subprocess
        import types

        import gside.analysis.taxonomic_validator as tvmod

        out = tmp_path / "o"
        self._checkm2_tsv(out, ["g1\t99.0\t1.0"])
        monkeypatch.setattr(tvmod, "CHECKM2_DB", tmp_path)
        monkeypatch.setattr(tvmod, "_find_tool", lambda name: "/bin/checkm2")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        assert tvmod._run_checkm2("q.fna", out) == (99.0, 1.0)

    def test_checkm2_bad_floats(self, tmp_path, monkeypatch):
        import subprocess
        import types

        import gside.analysis.taxonomic_validator as tvmod

        out = tmp_path / "o"
        self._checkm2_tsv(out, ["g1\tbad\talsobad"])
        monkeypatch.setattr(tvmod, "CHECKM2_DB", tmp_path)
        monkeypatch.setattr(tvmod, "_find_tool", lambda name: "/bin/checkm2")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        assert tvmod._run_checkm2("q.fna", out) == (None, None)

    def test_checkm2_short_file(self, tmp_path, monkeypatch):
        import subprocess
        import types

        import gside.analysis.taxonomic_validator as tvmod

        out = tmp_path / "o"
        out.mkdir(parents=True, exist_ok=True)
        (out / "checkm2_results.tsv").write_text("Name\tCompleteness\n")
        monkeypatch.setattr(tvmod, "CHECKM2_DB", tmp_path)
        monkeypatch.setattr(tvmod, "_find_tool", lambda name: "/bin/checkm2")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        assert tvmod._run_checkm2("q.fna", out) == (None, None)

    def test_gtdbtk_ar122(self, tmp_path, monkeypatch):
        import subprocess
        import types

        import gside.analysis.taxonomic_validator as tvmod

        out = tmp_path / "o"
        (out / "gtdb_output").mkdir(parents=True)
        (out / "gtdb_output" / "gtdbtk.ar122.summary.tsv").write_text(
            "user_genome\tclassification\nq\td__Archaea;s__X\n"
        )
        monkeypatch.setattr(tvmod, "GTDB_DB", tmp_path)
        monkeypatch.setattr(tvmod, "_find_tool", lambda name: "/bin/gtdbtk")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        q = tmp_path / "q.fna"
        q.write_text(">q\nACGT\n")
        assert tvmod._run_gtdbtk(str(q), out) == "d__Archaea;s__X"

    def test_gtdbtk_short_file(self, tmp_path, monkeypatch):
        import subprocess
        import types

        import gside.analysis.taxonomic_validator as tvmod

        out = tmp_path / "o"
        (out / "gtdb_output").mkdir(parents=True)
        (out / "gtdb_output" / "gtdbtk.bac120.summary.tsv").write_text("user_genome\n")
        monkeypatch.setattr(tvmod, "GTDB_DB", tmp_path)
        monkeypatch.setattr(tvmod, "_find_tool", lambda name: "/bin/gtdbtk")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        q = tmp_path / "q.fna"
        q.write_text(">q\nACGT\n")
        assert tvmod._run_gtdbtk(str(q), out) == ""
