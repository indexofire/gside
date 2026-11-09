"""Paired-end FASTQ auto-grouping: mates must form one sample, not two."""

from __future__ import annotations

import types

from click.testing import CliRunner

from gside import cli
from gside.analysis import read_screening as rs
from gside.cli import main


class TestIsolateDispatch:
    def test_pair_is_one_sample(self, monkeypatch):
        seen: list[list[str]] = []

        def fake_run_species(samples, *args, **kwargs):
            seen.extend(samples)
            return 0

        # _run_species receives grouped samples (patched at the dispatch edge)
        monkeypatch.setattr(cli, "_run_species", fake_run_species)
        result = CliRunner().invoke(main, ["scan", "S_R1.fq", "S_R2.fq", "-m", "marker"])
        assert result.exit_code == 0
        assert seen == [["S_R1.fq", "S_R2.fq"]]

    def test_two_pairs_are_two_samples(self, monkeypatch):
        seen: list[list[str]] = []
        monkeypatch.setattr(cli, "_run_species", lambda samples, *a, **k: seen.extend(samples) or 0)
        result = CliRunner().invoke(
            main, ["scan", "A_R1.fq", "A_R2.fq", "B_R1.fq", "B_R2.fq", "-m", "marker"]
        )
        assert result.exit_code == 0
        assert seen == [["A_R1.fq", "A_R2.fq"], ["B_R1.fq", "B_R2.fq"]]


class TestIdentifyOneFiles:
    def test_marker_receives_both_mates(self, monkeypatch):
        import gside.analysis.multigene_identifier as mg

        seen: dict[str, object] = {}

        def fake(files, threads=4):
            seen["files"] = files
            return types.SimpleNamespace(
                to_dict=lambda: {"species": "Salmonella", "confidence": "high", "database": {}}
            )

        monkeypatch.setattr(mg, "identify_multigene", fake)
        payload = cli._identify_one(["S_R1.fq", "S_R2.fq"], "marker", None)
        assert seen["files"] == ["S_R1.fq", "S_R2.fq"]
        assert payload["contigs"] == "S_R1.fq"
        assert payload["methods"]["marker"]["species"] == "Salmonella"

    def test_sourmash_receives_both_mates(self, monkeypatch):
        import gside.analysis.sourmash_identifier as sm

        seen: dict[str, object] = {}
        monkeypatch.setattr(
            sm,
            "identify_by_sourmash",
            lambda files, db_dir=None, dbs=None: (
                seen.update(files=files)
                or types.SimpleNamespace(
                    to_dict=lambda: {"species": "X", "confidence": "high", "database": {}}
                )
            ),
        )
        cli._identify_one(["S_R1.fq", "S_R2.fq"], "sourmash", None)
        assert seen["files"] == ["S_R1.fq", "S_R2.fq"]

    def test_mash_receives_both_mates(self, monkeypatch):
        import gside.analysis.ani_identifier as ani

        seen: dict[str, object] = {}
        monkeypatch.setattr(
            ani,
            "identify_by_ani",
            lambda files, mode, db_dir=None, threads=4, dbs=None: (
                seen.update(files=files)
                or types.SimpleNamespace(
                    to_dict=lambda: {"species": "X", "confidence": "high", "database": {}}
                )
            ),
        )
        cli._identify_one(["S_R1.fq", "S_R2.fq"], "mash_refseq", None)
        assert seen["files"] == ["S_R1.fq", "S_R2.fq"]


class TestScreenReadsMultipleFiles:
    def test_command_lists_both_files(self, monkeypatch):
        seen = {}

        def fake_run(cmd, **kwargs):
            seen["cmd"] = cmd
            return types.SimpleNamespace(stdout="")

        monkeypatch.setattr(rs, "require_bin", lambda *a, **k: "/bin/minimap2")
        monkeypatch.setattr(rs, "run_checked", fake_run)
        monkeypatch.setattr(rs, "mean_read_length", lambda path: 150.0)
        rs.screen_reads(["S_R1.fq", "S_R2.fq"], "markers.fasta")
        assert seen["cmd"][-3:] == ["markers.fasta", "S_R1.fq", "S_R2.fq"]
