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
        assert calls["args"][0] == ["g.fna"]
        assert calls["args"][1] == "marker"

    def test_fastq_marker_is_allowed(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["scan", "r.fq", "-m", "marker"]).exit_code == 0
        assert calls["args"][1] == "marker"

    def test_fastq_rejected_for_genome_tiers(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "panel"])
        assert result.exit_code == 2
        assert "does not accept FASTQ" in result.output
        assert "add --meta" not in result.output

    def test_fastq_hint_mentions_meta(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "panel"])
        assert "supported: FASTA" in result.output

    def test_meta_only_mode_requires_flag(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "-m", "kraken2"])
        assert result.exit_code == 2
        assert "add --meta" in result.output

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

    def test_samples_require_meta(self, tmp_path):
        manifest = tmp_path / "s.tsv"
        manifest.write_text("A\ta.fq\n")
        result = CliRunner().invoke(main, ["scan", "g.fna", "--samples", str(manifest)])
        assert result.exit_code == 2
        assert "--samples is only available with --meta" in result.output

    def test_batch_jobs_forwarded(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_species")
        assert CliRunner().invoke(main, ["scan", "a.fna", "b.fna", "-j", "4"]).exit_code == 0
        assert calls["args"][5] == 4


class TestScanMeta:
    def test_fasta_kraken2(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        assert CliRunner().invoke(main, ["scan", "g.fna", "--meta", "-m", "kraken2"]).exit_code == 0
        assert calls["args"][0] == ["g.fna"]
        assert calls["args"][1] == "kraken2"

    def test_fastq_sylph(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        assert CliRunner().invoke(main, ["scan", "r.fq", "--meta", "-m", "sylph"]).exit_code == 0
        assert calls["args"][1] == "sylph"

    def test_fasta_sylph_rejected(self):
        result = CliRunner().invoke(main, ["scan", "g.fna", "--meta", "-m", "sylph"])
        assert result.exit_code == 2
        assert "does not accept FASTA" in result.output

    def test_isolate_mode_with_meta_rejected(self):
        result = CliRunner().invoke(main, ["scan", "g.fna", "--meta", "-m", "marker"])
        assert result.exit_code == 2
        assert "not available with --meta" in result.output
        assert "drop --meta" in result.output

    def test_default_mode_is_all(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        assert CliRunner().invoke(main, ["scan", "r.fq", "--meta"]).exit_code == 0
        assert calls["args"][1] == "all"

    def test_fasta_all_rejected_because_sylph_needs_reads(self):
        result = CliRunner().invoke(main, ["scan", "g.fna", "--meta", "-m", "all"])
        assert result.exit_code == 2
        assert "sylph" in result.output and "FASTQ" in result.output
        assert "use -m kraken2" in result.output

    def test_fasta_kraken2_explicit(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        assert CliRunner().invoke(main, ["scan", "g.fna", "--meta", "-m", "kraken2"]).exit_code == 0
        assert calls["args"][1] == "kraken2"

    def test_paired_fastq(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        assert CliRunner().invoke(main, ["scan", "a.fq", "b.fq", "--meta"]).exit_code == 0
        assert calls["args"][0] == ["a.fq", "b.fq"]

    def test_paired_mixed_input_rejected(self):
        result = CliRunner().invoke(main, ["scan", "a.fq", "b.fna", "--meta"])
        assert result.exit_code == 2
        assert "paired input must be two FASTQ" in result.output

    def test_three_files_require_samples(self):
        result = CliRunner().invoke(main, ["scan", "a.fq", "b.fq", "c.fq", "--meta"])
        assert result.exit_code == 2
        assert "at most 2 files" in result.output

    def test_jobs_requires_samples(self):
        result = CliRunner().invoke(main, ["scan", "r.fq", "--meta", "-j", "2"])
        assert result.exit_code == 2
        assert "--jobs requires --samples" in result.output

    def test_manifest_dispatch(self, monkeypatch, tmp_path):
        calls = _capture(monkeypatch, "_run_meta_samples")
        manifest = tmp_path / "s.tsv"
        manifest.write_text("A\ta.fq\nB\tb.fq\n")
        result = CliRunner().invoke(main, ["scan", "--meta", "--samples", str(manifest)])
        assert result.exit_code == 0
        samples = calls["args"][0]
        assert [s["id"] for s in samples] == ["A", "B"]
        assert calls["args"][1] == "all"

    def test_manifest_validates_each_sample(self, monkeypatch, tmp_path):
        _capture(monkeypatch, "_run_meta_samples")
        manifest = tmp_path / "s.tsv"
        manifest.write_text("A\ta.fq\nB\tb.fna\n")
        result = CliRunner().invoke(
            main, ["scan", "--meta", "-m", "sylph", "--samples", str(manifest)]
        )
        assert result.exit_code == 2
        assert "does not accept FASTA" in result.output

    def test_files_and_samples_are_exclusive(self, tmp_path):
        manifest = tmp_path / "s.tsv"
        manifest.write_text("A\ta.fq\n")
        result = CliRunner().invoke(main, ["scan", "r.fq", "--meta", "--samples", str(manifest)])
        assert result.exit_code == 2
        assert "not both" in result.output


class TestAliases:
    def test_species_alias_forwards_without_meta(self, monkeypatch):
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
        assert "not available without --meta" in result.output

    def test_meta_alias_forces_meta(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        assert CliRunner().invoke(main, ["meta", "r.fq", "-m", "kraken2"]).exit_code == 0
        assert calls["args"][1] == "kraken2"

    def test_meta_alias_defaults_to_all(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        assert CliRunner().invoke(main, ["meta", "r.fq"]).exit_code == 0
        assert calls["args"][1] == "all"

    def test_meta_alias_rejects_isolate_modes(self):
        result = CliRunner().invoke(main, ["meta", "g.fna", "-m", "panel"])
        assert result.exit_code == 2
        assert "not available with --meta" in result.output

    def test_scan_is_the_only_visible_command(self):
        result = CliRunner().invoke(main, ["--help"])
        commands = result.output.split("Commands:", 1)[1]
        assert "  scan" in commands
        assert "  species" not in commands
        assert "  meta" not in commands

    def test_meta_alias_help_still_works(self):
        result = CliRunner().invoke(main, ["meta", "-h"])
        assert result.exit_code == 0
        assert "-m, --mode" in result.output


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

    def test_meta_combination_is_forwarded(self, monkeypatch):
        calls = _capture(monkeypatch, "_run_meta")
        result = CliRunner().invoke(main, ["scan", "r.fq", "--meta", "-m", "kraken2+sylph"])
        assert result.exit_code == 0
        assert calls["args"][1] == "kraken2+sylph"

    def test_meta_combination_rejects_fasta_component(self):
        result = CliRunner().invoke(main, ["scan", "g.fna", "--meta", "-m", "kraken2+sylph"])
        assert result.exit_code == 2
        assert "mode 'sylph' does not accept FASTA" in result.output


class TestExpandModes:
    def test_all_expands_per_context(self):
        assert cli._expand_modes("all", meta=False) == [
            "marker",
            "panel",
            "mash_refseq",
            "sourmash",
            "gtdbtk",
        ]
        assert cli._expand_modes("all", meta=True) == ["kraken2", "sylph"]

    def test_duplicates_are_dropped_in_order(self):
        assert cli._expand_modes("marker+panel+marker", meta=False) == ["marker", "panel"]

    def test_whitespace_is_tolerated(self):
        assert cli._expand_modes("marker + sourmash", meta=False) == ["marker", "sourmash"]
