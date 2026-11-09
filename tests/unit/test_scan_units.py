"""Unified `gside scan` command: mode/input validation and alias behaviour."""

from __future__ import annotations

from typing import Any

import pytest
from click.testing import CliRunner

from gside import cli
from gside.cli import _input_kind, main


def _capture(monkeypatch, name: str) -> dict[str, Any]:
    calls: dict[str, Any] = {}

    def fake(*args: Any, **kwargs: Any) -> int:
        calls["args"] = args
        calls["kwargs"] = kwargs
        return 0

    monkeypatch.setattr(cli, name, fake)
    return calls


class TestInputKind:
    @pytest.mark.parametrize(
        "name,kind",
        [
            ("g.fna", "fasta"),
            ("g.fa", "fasta"),
            ("g.fasta", "fasta"),
            ("g.fna.gz", "fasta"),
            ("r.fq", "fastq"),
            ("r.fastq", "fastq"),
            ("r.fq.gz", "fastq"),
            ("r.FASTQ.GZ", "fastq"),
            ("notes.txt", "unknown"),
            ("noext", "unknown"),
        ],
    )
    def test_extension_matrix(self, name, kind):
        assert _input_kind(name) == kind


class TestScanIsolate:
    def test_defaults_to_marker(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["scan", "g.fna"]).exit_code == 0
        assert calls["args"][0] == [["g.fna"]]
        assert calls["args"][1] == "marker"

    def test_fastq_marker_is_allowed(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["scan", "r.fq", "-m", "marker"]).exit_code == 0
        assert calls["args"][1] == "marker"

    def test_fastq_rejected_for_genome_tiers(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "panel"])
        assert result.exit_code == 2
        assert "does not accept FASTQ" in result.output

    def test_all_on_fastq_is_rejected(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "all"])
        assert result.exit_code == 2
        assert "does not accept FASTQ" in result.output

    def test_unknown_extension_is_rejected(self):
        result = CliRunner().invoke(main, ["scan", "notes.txt", "-m", "marker"])
        assert result.exit_code == 2
        assert "cannot infer input type" in result.output

    def test_missing_input(self):
        result = CliRunner().invoke(main, ["scan"])
        assert result.exit_code == 2
        assert "missing input" in result.output

    def test_batch_jobs_forwarded(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["scan", "a.fna", "b.fna", "-j", "4"]).exit_code == 0
        assert calls["args"][5] == 4


class TestMetaMovedToGmeta:
    """The meta layer (kraken2/sylph community profiling) lives in gmeta now."""

    def test_meta_flag_is_gone(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "--meta"])
        assert result.exit_code == 2

    def test_kraken2_mode_is_unknown(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "kraken2"])
        assert result.exit_code == 2
        assert "unknown mode" in result.output

    def test_sylph_mode_is_unknown(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "sylph"])
        assert result.exit_code == 2
        assert "unknown mode" in result.output

    def test_samples_option_is_gone(self, tmp_path):
        manifest = tmp_path / "s.tsv"
        manifest.write_text("A\ta.fq\n")
        result = CliRunner().invoke(main, ["scan", "g.fna", "--samples", str(manifest)])
        assert result.exit_code == 2

    def test_meta_alias_is_gone(self):
        result = CliRunner().invoke(main, ["meta", "r.fq"])
        assert result.exit_code != 0


class TestAliases:
    def test_species_alias_forwards(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["species", "g.fna", "-m", "panel"]).exit_code == 0
        assert calls["args"][1] == "panel"

    def test_species_alias_defaults_to_marker(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["species", "g.fna"]).exit_code == 0
        assert calls["args"][1] == "marker"

    def test_species_alias_rejects_meta_modes(self):
        result = CliRunner().invoke(main, ["species", "r.fq", "-m", "kraken2"])
        assert result.exit_code == 2
        assert "unknown mode" in result.output

    def test_scan_is_the_only_visible_command(self):
        result = CliRunner().invoke(main, ["--help"])
        commands = result.output.split("Commands:", 1)[1]
        assert "  scan" in commands
        assert "  species" not in commands
        assert "  meta" not in commands


class TestReadsCapableGenomeTiers:
    def test_mash_refseq_accepts_fastq(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["scan", "r.fq", "-m", "mash_refseq"]).exit_code == 0
        assert calls["args"][1] == "mash_refseq"

    def test_sourmash_accepts_fastq(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["scan", "r.fq", "-m", "sourmash"]).exit_code == 0
        assert calls["args"][1] == "sourmash"

    def test_panel_still_rejects_fastq(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "panel"])
        assert result.exit_code == 2
        assert "does not accept FASTQ" in result.output

    def test_all_still_rejects_fastq(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "all"])
        assert result.exit_code == 2
        assert "does not accept FASTQ" in result.output


class TestCompositeModes:
    def test_isolate_combination_is_forwarded(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        result = CliRunner().invoke(main, ["scan", "g.fna", "-m", "marker+sourmash"])
        assert result.exit_code == 0
        assert calls["args"][1] == "marker+sourmash"

    def test_combination_of_reads_capable_tiers(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "marker+mash_refseq+sourmash"])
        assert result.exit_code == 0
        assert calls["args"][1] == "marker+mash_refseq+sourmash"

    def test_incompatible_component_is_named(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "marker+panel"])
        assert result.exit_code == 2
        assert "mode 'panel' does not accept FASTQ" in result.output

    def test_all_cannot_be_combined(self):
        result = CliRunner().invoke(main, ["scan", "g.fna", "-m", "all+marker"])
        assert result.exit_code == 2
        assert "cannot be combined" in result.output

    def test_empty_mode_rejected(self):
        result = CliRunner().invoke(main, ["scan", "g.fna", "-m", "  "])
        assert result.exit_code == 2
        assert "--mode is empty" in result.output


class TestExpandModes:
    def test_all_expands_to_isolate_modes(self):
        assert cli._expand_modes("all") == [
            "marker",
            "panel",
            "mash_refseq",
            "sourmash",
            "gtdbtk",
        ]

    def test_duplicates_are_dropped_in_order(self):
        assert cli._expand_modes("marker+panel+marker") == ["marker", "panel"]

    def test_whitespace_is_tolerated(self):
        assert cli._expand_modes("marker + sourmash") == ["marker", "sourmash"]


def _marker_family_db(root, name="food"):
    path = root / "marker" / name
    path.mkdir(parents=True)
    (path / "markers.fasta").write_text(">markers~~~g1~~~x marker\nACGTACGT\n")
    (path / "marker_rules.yaml").write_text("rules:\n- species: X\n  genes: [g1]\n")
    return path


def _panel_family_db(root, name="lab2026"):
    path = root / "panel" / name
    (path / "panel.sketch").mkdir(parents=True)
    (path / "panel.sketch" / "sketches.db").write_bytes(b"sketch")
    (path / "metadata.tsv").write_text("genome\tspecies\n")
    return path


class TestDbSelection:
    def test_named_selection_reaches_run_species(self, tmp_path, monkeypatch):
        db = _marker_family_db(tmp_path)
        captured: dict[str, Any] = {}
        monkeypatch.setattr(cli, "_run_species", lambda *a, **k: captured.update(kwargs=k) or 0)
        result = CliRunner().invoke(main, ["scan", "g.fna", "-D", str(tmp_path), "-d", "food"])
        assert result.exit_code == 0
        assert captured["kwargs"]["db_selections"] == {"marker": [db]}

    def test_marker_db_alias_folds_into_selections(self, tmp_path):
        db = _marker_family_db(tmp_path)
        captured: dict[str, Any] = {}
        orig = cli._run_species
        cli._run_species = lambda *a, **k: captured.update(kwargs=k) or 0
        try:
            result = CliRunner().invoke(
                main, ["scan", "g.fna", "-D", str(tmp_path), "--marker-db", str(db)]
            )
        finally:
            cli._run_species = orig
        assert result.exit_code == 0
        assert captured["kwargs"]["db_selections"] == {"marker": [db]}

    def test_mode_qualified_multi_mode(self, tmp_path):
        marker = _marker_family_db(tmp_path)
        panel = _panel_family_db(tmp_path)
        captured: dict[str, Any] = {}
        orig = cli._run_species
        cli._run_species = lambda *a, **k: captured.update(kwargs=k) or 0
        try:
            result = CliRunner().invoke(
                main,
                [
                    "scan",
                    "g.fna",
                    "-m",
                    "marker+panel",
                    "-D",
                    str(tmp_path),
                    "-d",
                    "marker:food",
                    "-d",
                    "panel:lab2026",
                ],
            )
        finally:
            cli._run_species = orig
        assert result.exit_code == 0
        assert captured["kwargs"]["db_selections"] == {"marker": [marker], "panel": [panel]}

    def test_unknown_mode_prefix_is_usage_error(self, tmp_path):
        _marker_family_db(tmp_path)
        result = CliRunner().invoke(
            main, ["scan", "g.fna", "-D", str(tmp_path), "-d", "panel:lab2026"]
        )
        assert result.exit_code == 2
        assert "not part of this invocation" in result.output

    def test_bare_name_with_multiple_modes_is_usage_error(self, tmp_path):
        _marker_family_db(tmp_path)
        result = CliRunner().invoke(
            main,
            ["scan", "g.fna", "-m", "marker+panel", "-D", str(tmp_path), "-d", "food"],
        )
        assert result.exit_code == 2
        assert "qualify the database" in result.output

    def test_mode_without_family_is_usage_error(self, tmp_path):
        result = CliRunner().invoke(main, ["scan", "g.fna", "-m", "gtdbtk", "-d", "whatever"])
        assert result.exit_code == 2
        assert "no mode in this invocation supports named databases" in result.output
