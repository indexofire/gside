"""Pure unit tests for cli.py (no binaries, no databases, milliseconds)."""

from __future__ import annotations

import json
import logging

import pytest
from click.testing import CliRunner

from gside import __version__, cli
from gside.cli import (
    TABLE_COLUMNS,
    _identify_one,
    _method_cell,
    _render_table,
    _run_species,
    main,
    setup_logging,
)


def _payload(contigs="a.fna", species="Salmonella", conf="high", rule="inva (1/1 required)"):
    return {
        "contigs": contigs,
        "methods": {
            "marker": {"species": species, "confidence": conf, "matched_rule": rule},
        },
        "verdict": {"species": species, "confidence": conf, "basis": ["marker"]},
    }


def _method_payload(contigs="a.fna", method="sourmash", species="S", conf="high"):
    return {
        "contigs": contigs,
        "methods": {method: {"result": {"species": species, "confidence": conf}}},
        "verdict": {"species": species, "confidence": conf, "basis": [method]},
    }


def _missing_db_payload(contigs="q.fna", method="sourmash", tier=None, size="~3.9GB"):
    tier = method if tier is None else tier
    return {
        "contigs": contigs,
        "methods": {
            method: {
                "method": method,
                "error": f"{tier} database not found in /x (run: gside db setup --tier {tier})",
                "missing_tier": tier,
                "suggested_command": f"gside db setup --tier {tier}",
                "estimated_size": size,
            }
        },
        "verdict": {"species": "Unknown", "confidence": "low", "basis": [method]},
    }


class TestSetupLogging:
    def test_conflict_raises(self):
        with pytest.raises(ValueError, match="verbose and quiet"):
            setup_logging(verbose=True, quiet=True)

    def test_levels(self, monkeypatch):
        calls = {}
        monkeypatch.setattr(logging, "basicConfig", lambda **kw: calls.update(kw))
        setup_logging(verbose=True)
        assert calls["level"] == logging.DEBUG
        setup_logging(quiet=True)
        assert calls["level"] == logging.ERROR
        setup_logging()
        assert calls["level"] == logging.WARNING


class TestMethodCell:
    def test_missing_mode(self):
        assert _method_cell({}, "panel", "species") == ""

    def test_nested_result(self):
        methods = {"panel": {"result": {"species": "X", "confidence": "high"}}}
        assert _method_cell(methods, "panel", "species") == "X"

    def test_none_becomes_empty(self):
        assert _method_cell({"m": {"species": None}}, "m", "species") == ""


class TestRenderTable:
    def test_multi_row_and_sanitize(self):
        evil = _payload(contigs="a|b\nc.fna")
        rows = _render_table([self._payload_row(), evil], markdown=False).splitlines()
        assert rows[0].split("\t") == TABLE_COLUMNS
        assert len(rows) == 3
        assert "|" not in rows[2] and "\n" not in rows[2]

    def test_md_column_count(self):
        out = _render_table([self._payload_row()], markdown=True).splitlines()
        assert out[1].count("---") == len(TABLE_COLUMNS)

    def test_tab_and_cr_in_cells_do_not_add_columns(self):
        evil = _payload(contigs="a\tb\rc.fna")
        row = _render_table([evil], markdown=False).splitlines()[1].split("\t")
        assert len(row) == len(TABLE_COLUMNS)
        assert "\r" not in row[0]

    def test_empty_methods(self):
        payload = {"contigs": "x.fna", "methods": {}, "verdict": {}}
        row = _render_table([payload], markdown=False).splitlines()[1].split("\t")
        assert len(row) == len(TABLE_COLUMNS)
        assert row[1] == ""

    @staticmethod
    def _payload_row():
        return _payload()


class TestRunSpeciesBatch:
    def test_json_single_is_object(self, monkeypatch, capsys):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        assert _run_species(["a.fna"], "marker", None) == 0
        assert json.loads(capsys.readouterr().out)["contigs"] == "a.fna"

    def test_json_multi_is_array(self, monkeypatch, capsys):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        assert _run_species(["a.fna", "b.fna"], "marker", None) == 0
        docs = json.loads(capsys.readouterr().out)
        assert [d["contigs"] for d in docs] == ["a.fna", "b.fna"]

    def test_tsv_and_md(self, monkeypatch, capsys):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        assert _run_species(["a.fna"], "marker", None, "tsv") == 0
        assert capsys.readouterr().out.startswith("contigs\t")
        assert _run_species(["a.fna"], "marker", None, "md") == 0
        assert capsys.readouterr().out.startswith("| contigs |")

    def test_identify_one_missing_file_surfaces_error(self, tmp_path, capsys):
        missing = str(tmp_path / "nope.fna")
        assert _run_species([missing], "marker", None) == 0
        payload = json.loads(capsys.readouterr().out)
        assert "error" in payload["methods"]["marker"]
        assert payload["verdict"]["species"] == "Unknown"
        assert payload["verdict"]["confidence"] == "low"

    def test_identify_one_mash_missing_db_surfaces_error(self, tmp_path):
        payload = _identify_one("q.fna", "mash_refseq", db_dir=str(tmp_path))
        assert "gside db setup --tier mash" in payload["methods"]["mash_refseq"]["error"]
        assert payload["verdict"]["species"] == "Unknown"
        assert payload["verdict"]["confidence"] == "low"


class TestMainGroup:
    def test_verbose_quiet_conflict(self):
        result = CliRunner().invoke(main, ["--verbose", "--quiet"])
        assert result.exit_code == 2
        assert "cannot be used together" in result.output

    def test_verbose_quiet_short_flags_conflict(self):
        result = CliRunner().invoke(main, ["-v", "-q"])
        assert result.exit_code == 2
        assert "cannot be used together" in result.output

    def test_bare_help_exits_2(self):
        result = CliRunner().invoke(main, [])
        assert result.exit_code == 2
        assert "Usage" in result.output

    def test_verbose_only_no_subcommand_shows_help(self):
        result = CliRunner().invoke(main, ["-v"])
        assert result.exit_code == 0
        assert "Usage" in result.output


class TestMultigeneBranches:
    def test_near_threshold_note(self, monkeypatch):
        import gside.analysis.multigene_identifier as mg

        monkeypatch.setattr(
            mg,
            "_blast_contigs",
            lambda contigs, threads=4: [
                {"gene": "inva", "identity": 80.0, "coverage": 50.0},
            ],
        )
        r = mg.identify_multigene("dummy.fna").to_dict()
        assert r["species"] == "Unknown"
        assert any("Near-threshold" in n for n in r["notes"])

    def test_no_hits_no_notes(self, monkeypatch):
        import gside.analysis.multigene_identifier as mg

        monkeypatch.setattr(mg, "_blast_contigs", lambda contigs, threads=4: [])
        r = mg.identify_multigene("dummy.fna").to_dict()
        assert r["species"] == "Unknown"
        assert r["notes"] == []


class TestDbCommands:
    def _db(self, tmp_path, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "SPECIES_DB_DIR", tmp_path / "db")
        return dbmod

    def test_status(self, tmp_path, monkeypatch, capsys):
        from click.testing import CliRunner

        self._db(tmp_path, monkeypatch)
        result = CliRunner().invoke(main, ["db", "status"])
        assert result.exit_code == 0
        assert "D1_marker" in result.output

    def test_list(self, tmp_path, monkeypatch, capsys):
        from click.testing import CliRunner

        self._db(tmp_path, monkeypatch)
        result = CliRunner().invoke(main, ["db", "list"])
        assert result.exit_code == 0
        assert "sourmash" in result.output

    def test_bare_db_invokes_status(self, monkeypatch):
        from click.testing import CliRunner

        import gside.db as dbmod

        seen = []
        monkeypatch.setattr(dbmod, "run_db_command", lambda args: seen.append(args) or 0)
        result = CliRunner().invoke(main, ["db"])
        assert result.exit_code == 0
        assert seen == [["status"]]

    def test_setup_ok_and_error(self, tmp_path, monkeypatch):
        from click.testing import CliRunner

        import gside.db as dbmod

        self._db(tmp_path, monkeypatch)
        (tmp_path / "db" / "D2_ani" / "panel.sketch").mkdir(parents=True)
        (tmp_path / "db" / "D2_ani" / "panel.sketch" / "sketches.db").write_bytes(b"x")
        assert CliRunner().invoke(main, ["db", "setup", "--tier", "panel"]).exit_code == 0
        monkeypatch.setattr(dbmod, "db_setup", lambda tier, source: {"x": "ERROR: boom"})
        assert CliRunner().invoke(main, ["db", "setup", "--tier", "panel"]).exit_code == 1


class TestIdentifyModes:
    def test_all_branches(self, monkeypatch):
        import gside.analysis.ani_identifier as ani
        import gside.analysis.gtdbtk_identifier as gmod
        import gside.analysis.multigene_identifier as mg
        import gside.analysis.sourmash_identifier as smod
        import gside.cli as climod

        monkeypatch.setattr(mg, "identify_multigene", lambda c, threads=4: _FakeResult("M", "high"))
        monkeypatch.setattr(
            ani, "identify_by_ani", lambda c, mode, db_dir=None, threads=4: _FakeResult("A", "high")
        )
        monkeypatch.setattr(
            smod, "identify_by_sourmash", lambda c, db_dir=None: _FakeResult("S", "high")
        )
        monkeypatch.setattr(
            gmod, "identify_by_gtdbtk", lambda c, threads=4: _FakeResult("G", "high")
        )
        payload = climod._identify_one("x.fna", "all", None)
        assert set(payload["methods"]) == {"marker", "panel", "mash_refseq", "sourmash", "gtdbtk"}
        assert payload["verdict"]["basis"] == ["gtdbtk"]

    def test_branch_error_entry(self, monkeypatch):
        import gside.analysis.multigene_identifier as mg
        import gside.cli as climod

        def boom(c, threads=4):
            raise OSError("blast exploded")

        monkeypatch.setattr(mg, "identify_multigene", boom)
        payload = climod._identify_one("x.fna", "marker", None)
        assert payload["methods"]["marker"]["error"] == "blast exploded"
        assert payload["verdict"]["species"] == "Unknown"


class _FakeResult:
    def __init__(self, species, confidence):
        self._d = {"species": species, "confidence": confidence}

    def to_dict(self):
        return dict(self._d)


class TestCommandsEndToEnd:
    def test_species_command(self, monkeypatch):
        from click.testing import CliRunner

        import gside.cli as climod

        monkeypatch.setattr(climod, "_run_species", lambda *a: 0)
        assert CliRunner().invoke(main, ["species", "x.fna"]).exit_code == 0

    def test_validate_command(self, monkeypatch, tmp_path):
        from click.testing import CliRunner

        import gside.analysis.taxonomic_validator as tv

        seen: dict[str, object] = {}

        def fake(contigs, output_dir, threads=4):
            seen.update(contigs=contigs, output_dir=output_dir, threads=threads)
            return tv.ValidateResult(status="ok", quality={"verdict": "PASS"})

        monkeypatch.setattr(tv, "validate_genome", fake)
        source = tmp_path / "x.fna"
        source.write_text(">c\nACGT\n")
        result = CliRunner().invoke(
            main, ["validate", str(source), "-o", str(tmp_path / "tax"), "-t", "8"]
        )
        assert result.exit_code == 0, result.output
        assert seen["threads"] == 8
        assert str(seen["output_dir"]).endswith("tax")
        assert '"verdict": "PASS"' in result.output


class TestSpeciesShortFlags:
    def test_short_flags_parse_like_long(self, monkeypatch):
        seen = []
        monkeypatch.setattr(cli, "_run_species", lambda *a: seen.append(a) or 0)
        short = CliRunner().invoke(
            main, ["species", "-m", "marker", "-f", "tsv", "-D", "/db", "x.fna"]
        )
        long = CliRunner().invoke(
            main,
            ["species", "--mode", "marker", "--format", "tsv", "--db-dir", "/db", "x.fna"],
        )
        assert short.exit_code == 0
        assert long.exit_code == 0
        assert seen[0] == seen[1]


class TestSpeciesOutputOption:
    def test_output_writes_file_stdout_clean_stderr_notice(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        out = tmp_path / "out.json"
        result = CliRunner().invoke(main, ["species", "--output", str(out), "a.fna"])
        assert result.exit_code == 0
        assert json.loads(out.read_text())["contigs"] == "a.fna"
        assert result.stdout == ""
        assert f"Results written to {out}" in result.stderr

    def test_output_short_flag_multi_file_array(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        out = tmp_path / "o.json"
        result = CliRunner().invoke(main, ["species", "-o", str(out), "a.fna", "b.fna"])
        assert result.exit_code == 0
        assert [d["contigs"] for d in json.loads(out.read_text())] == ["a.fna", "b.fna"]
        assert "Results written to" in result.stderr

    def test_output_tsv_payload(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        out = tmp_path / "o.tsv"
        result = CliRunner().invoke(main, ["species", "-f", "tsv", "-o", str(out), "a.fna"])
        assert result.exit_code == 0
        assert out.read_text().startswith("contigs\t")
        assert result.stdout == ""


class TestDbOutputFormats:
    _STATUS = {
        "D1_marker": {"ready": True, "path": "/p/D1_marker", "tier": "mini"},
        "D2_ani": {"ready": False, "path": "/p/D2_ani", "tier": "panel"},
    }

    _TIERS = {
        "mini": {
            "description": "Marker rules + sequences (bundled)",
            "size": "~1MB",
            "items": ["D1_marker"],
        },
    }

    def test_status_tsv_machine_readable(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "status", "--format", "tsv"])
        assert result.exit_code == 0
        lines = result.output.splitlines()
        assert lines[0] == "COMPONENT\tREADY\tTIER\tPATH"
        assert lines[1] == "D1_marker\tready\tmini\t/p/D1_marker"
        assert lines[2] == "D2_ani\tmissing\tpanel\t/p/D2_ani"
        assert "✅" not in result.output
        assert "db setup" not in result.output

    def test_status_tsv_short_flag(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "status", "-f", "tsv"])
        assert result.exit_code == 0
        assert result.output.splitlines()[0] == "COMPONENT\tREADY\tTIER\tPATH"

    def test_status_json_roundtrip(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "status", "--format", "json"])
        assert result.exit_code == 0
        assert json.loads(result.output) == self._STATUS

    def test_status_table_default_has_borders_and_hint(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "status"])
        assert result.exit_code == 0
        assert "gside database status" in result.output
        assert "│ COMPONENT" in result.output
        assert "┌" in result.output and "└" in result.output
        assert "ready" in result.output and "missing" in result.output
        assert "gside db setup --tier <tier>" in result.output

    def test_list_tsv(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "_TIERS", self._TIERS)
        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "list", "--format", "tsv"])
        assert result.exit_code == 0
        assert result.output.splitlines() == [
            "NAME\tREADY\tSIZE\tDESCRIPTION",
            "mini\tready\t~1MB\tMarker rules + sequences (bundled)",
        ]

    def test_list_json(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "_TIERS", self._TIERS)
        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "list", "-f", "json"])
        assert result.exit_code == 0
        assert json.loads(result.output) == {"mini": {**self._TIERS["mini"], "ready": True}}


class TestDbSetupShortFlags:
    def test_short_tier_and_source(self, tmp_path, monkeypatch):
        import gside.db as dbmod

        seen = {}
        monkeypatch.setattr(dbmod, "SPECIES_DB_DIR", tmp_path / "db")
        monkeypatch.setattr(
            dbmod, "db_setup", lambda tier, source="": seen.update(tier=tier, source=source) or {}
        )
        result = CliRunner().invoke(main, ["db", "setup", "-t", "mash", "-s", str(tmp_path)])
        assert result.exit_code == 0
        assert seen == {"tier": "mash", "source": str(tmp_path)}


class TestHelpSurface:
    def test_main_help_lists_global_flags(self):
        result = CliRunner().invoke(main, ["--help"])
        assert result.exit_code == 0
        for token in ("-v, --verbose", "-q, --quiet", "-V, --version"):
            assert token in result.output

    def test_species_help_lists_short_flags(self):
        result = CliRunner().invoke(main, ["species", "-h"])
        assert result.exit_code == 0
        for token in ("-m, --mode", "-f, --format", "-D, --db-dir", "-o, --output"):
            assert token in result.output

    def test_db_status_and_list_help_show_format(self):
        for cmd in ("status", "list"):
            result = CliRunner().invoke(main, ["db", cmd, "-h"])
            assert result.exit_code == 0
            assert "-f, --format" in result.output
            assert "[default: table]" in result.output

    def test_db_setup_help_lists_short_flags(self):
        result = CliRunner().invoke(main, ["db", "setup", "-h"])
        assert result.exit_code == 0
        assert "-t, --tier" in result.output
        assert "-s, --source" in result.output

    def test_validate_inherits_h(self):
        result = CliRunner().invoke(main, ["validate", "-h"])
        assert result.exit_code == 0
        assert "-h, --help" in result.output


class TestRemainingBranches:
    def test_arbitrate_all_errors(self):
        from gside.cli import _arbitrate

        v = _arbitrate({"a": {"error": "x"}, "b": {"error": "y"}})
        assert v == {
            "species": "Unknown",
            "display_name": "Unknown",
            "strain": None,
            "confidence": "low",
            "basis": ["a", "b"],
            "supporting": [],
            "agreement": "0/0",
        }

    def test_setup_with_source(self, tmp_path, monkeypatch):
        from click.testing import CliRunner

        import gside.db as dbmod
        from gside.cli import main

        monkeypatch.setattr(dbmod, "SPECIES_DB_DIR", tmp_path / "db")
        result = CliRunner().invoke(
            main, ["db", "setup", "--tier", "panel", "--source", str(tmp_path)]
        )
        assert result.exit_code in (0, 1)


class TestVersionJson:
    def test_version_is_a_root_flag_not_a_subcommand(self):
        result = CliRunner().invoke(main, ["version"])
        assert result.exit_code == 2

    def test_version_text_unchanged(self):
        result = CliRunner().invoke(main, ["-V"])
        assert result.exit_code == 0
        assert result.output.strip() == f"gside, version {__version__}"

    def test_version_json_flag(self):
        result = CliRunner().invoke(main, ["-V", "-J"])
        assert result.exit_code == 0
        assert json.loads(result.output) == {"tool": "gside", "version": __version__}


class TestDbSetupJsonFormat:
    def test_json_mode_stdout_purity(self, monkeypatch):
        import gside.db as dbmod

        def fake_setup(tier, source=""):
            print("  downloading fake tier...")
            return {"D2_ani": "downloaded (fake)"}

        monkeypatch.setattr(dbmod, "db_setup", fake_setup)
        result = CliRunner().invoke(main, ["db", "setup", "-t", "panel", "-f", "json"])
        assert result.exit_code == 0
        assert json.loads(result.stdout) == {"D2_ani": "downloaded (fake)"}
        assert "downloading fake tier" in result.stderr
        assert "downloading fake tier" not in result.stdout

    def test_json_mode_error_still_exits_one(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "db_setup", lambda tier, source="": {"D2_ani": "ERROR: boom"})
        result = CliRunner().invoke(main, ["db", "setup", "-t", "panel", "-f", "json"])
        assert result.exit_code == 1
        assert json.loads(result.stdout) == {"D2_ani": "ERROR: boom"}


class TestStreamingOutput:
    def test_output_tsv_rows_appear_as_files_complete(self, monkeypatch, tmp_path):
        seen_during_run = {}

        def fake_identify(contigs, mode, db_dir, threads=4):
            if contigs == "b.fna":
                seen_during_run["mid"] = (tmp_path / "out.tsv").read_text()
            return _payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        cli._run_species(["a.fna", "b.fna"], "marker", None, "tsv", str(tmp_path / "out.tsv"))
        assert "a.fna" in seen_during_run["mid"].splitlines()[-1]
        assert "\t".join(TABLE_COLUMNS) in seen_during_run["mid"]

    def test_output_md_rows_stream(self, monkeypatch, tmp_path):
        seen_during_run = {}

        def fake_identify(contigs, mode, db_dir, threads=4):
            if contigs == "b.fna":
                seen_during_run["mid"] = (tmp_path / "out.md").read_text()
            return _payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        cli._run_species(["a.fna", "b.fna"], "marker", None, "md", str(tmp_path / "out.md"))
        assert "a.fna" in seen_during_run["mid"]
        assert "---" in seen_during_run["mid"]

    def test_stdout_rows_stream(self, monkeypatch, capsys):
        seen = {}

        def fake_identify(contigs, mode, db_dir, threads=4):
            if contigs == "b.fna":
                seen["mid"] = capsys.readouterr().out
            return _payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        cli._run_species(["a.fna", "b.fna"], "marker", None, "tsv", None)
        assert "a.fna" in seen["mid"]
        assert "b.fna" not in seen["mid"]

    def test_midbatch_failure_persists_prefix(self, monkeypatch, tmp_path):
        def fake_identify(contigs, mode, db_dir, threads=4):
            if contigs == "b.fna":
                raise RuntimeError("boom mid-batch")
            return _payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        with pytest.raises(RuntimeError):
            cli._run_species(["a.fna", "b.fna"], "marker", None, "tsv", str(tmp_path / "out.tsv"))
        persisted = (tmp_path / "out.tsv").read_text().splitlines()
        assert persisted[0].startswith("contigs")
        assert "a.fna" in persisted[1]
        assert len(persisted) == 2


class TestSpeciesParallelism:
    def test_jobs_and_threads_options_parse(self, monkeypatch):
        seen = {}

        def fake_run(*a, **k):
            seen.update(args=a, kwargs=k)
            return 0

        monkeypatch.setattr(cli, "_run_species", fake_run)
        result = CliRunner().invoke(main, ["species", "-j", "3", "-t", "8", "x.fna"])
        assert result.exit_code == 0
        assert seen["args"][5:] == (3, 8, False, False)

    def test_parallel_run_preserves_input_order(self, monkeypatch, tmp_path):
        import time

        done = []

        def fake_identify(contigs, mode, db_dir, threads=4):
            if contigs == "a.fna":
                time.sleep(0.25)
            done.append(contigs)
            return _payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        out = tmp_path / "o.tsv"
        cli._run_species(["a.fna", "b.fna", "c.fna"], "marker", None, "tsv", str(out), jobs=3)
        rows = out.read_text().splitlines()
        assert [r.split("\t")[0] for r in rows[1:]] == ["a.fna", "b.fna", "c.fna"]
        assert done[0] != "a.fna"

    def test_parallel_stream_holds_rows_until_input_turn(self, monkeypatch, tmp_path):
        import time

        seen = {}

        def fake_identify(contigs, mode, db_dir, threads=4):
            if contigs == "a.fna":
                time.sleep(0.25)
            if contigs == "c.fna":
                seen["mid"] = (tmp_path / "o.tsv").read_text()
            return _payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        cli._run_species(
            ["a.fna", "b.fna", "c.fna"], "marker", None, "tsv", str(tmp_path / "o.tsv"), jobs=3
        )
        assert "a.fna" not in seen["mid"]
        assert "b.fna" not in seen["mid"]

    def test_threads_flow_to_identifiers(self, monkeypatch):
        cap = {}

        class _Stub:
            def to_dict(self):
                return {}

        import gside.analysis.ani_identifier as ani
        import gside.analysis.multigene_identifier as mg

        monkeypatch.setattr(
            mg, "identify_multigene", lambda c, threads=4: cap.setdefault("mg", threads) or _Stub()
        )
        monkeypatch.setattr(
            ani,
            "identify_by_ani",
            lambda c, mode, db_dir=None, threads=4: cap.setdefault("ani", threads) or _Stub(),
        )
        cli._identify_one("q.fna", "marker", None, threads=8)
        cli._identify_one("q.fna", "panel", None, threads=8)
        assert cap["mg"] == 8 and cap["ani"] == 8


class TestSpeciesProgress:
    def test_progress_bar_module_exists(self):
        from gside.progress import IdleBatch, make_progress, species_progress

        assert callable(make_progress)
        assert callable(species_progress)
        idle = IdleBatch()
        idle.describe("x")
        idle.advance()

    def test_multi_file_not_quiet_not_tty_yields_idle(self, monkeypatch, tmp_path):
        from gside.progress import STATUS_CONSOLE, species_progress

        monkeypatch.setattr(type(STATUS_CONSOLE), "is_terminal", property(lambda self: False))
        with species_progress(["a", "b"], quiet=False) as (bar, tid):
            assert tid is None
            assert bar.advance() is None

    def test_quiet_yields_idle(self, monkeypatch):
        from gside.progress import species_progress

        with species_progress(["a", "b"], quiet=True) as (bar, tid):
            assert tid is None
            assert bar.describe("x") is None

    def test_single_file_yields_idle(self, monkeypatch):
        from gside.progress import species_progress

        with species_progress(["a"], quiet=False) as (bar, tid):
            assert tid is None

    def test_multi_file_tty_yields_live(self, monkeypatch):
        from gside.progress import STATUS_CONSOLE, species_progress

        monkeypatch.setattr(type(STATUS_CONSOLE), "is_terminal", property(lambda self: True))
        with species_progress(["a", "b", "c"], quiet=False) as (bar, tid):
            assert tid is not None
            bar.advance()


class TestNoHeaderOption:
    def test_no_header_flag_parses(self, monkeypatch):
        seen = {}
        monkeypatch.setattr(cli, "_run_species", lambda *a: seen.update(args=a) or 0)
        result = CliRunner().invoke(main, ["species", "--no-header", "x.fna"])
        assert result.exit_code == 0
        assert seen["args"][7] is True

    def test_tsv_no_header_omits_header_row(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        cli._run_species(["a.fna"], "marker", None, "tsv", None, no_header=True)
        # stdout was consumed by capsys in _run_species; verify via file
        cli._run_species(["a.fna"], "marker", None, "tsv", str(tmp_path / "o.tsv"), no_header=True)
        lines = (tmp_path / "o.tsv").read_text().strip().splitlines()
        assert not lines[0].startswith("contigs")
        assert "a.fna" in lines[0]

    def test_tsv_with_header_default(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        cli._run_species(["a.fna"], "marker", None, "tsv", str(tmp_path / "o.tsv"))
        lines = (tmp_path / "o.tsv").read_text().strip().splitlines()
        assert lines[0].split("\t")[0] == "contigs"

    def test_md_no_header_still_has_header(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        cli._run_species(["a.fna"], "marker", None, "md", str(tmp_path / "o.md"), no_header=True)
        lines = (tmp_path / "o.md").read_text().strip().splitlines()
        assert lines[0].startswith("| contigs")


class TestStructuredErrorFields:
    def test_sourmash_missing_db_gets_structured_fields(self, tmp_path):
        payload = _identify_one("q.fna", "sourmash", db_dir=str(tmp_path))
        method = payload["methods"]["sourmash"]
        assert method["missing_tier"] == "sourmash"
        assert method["suggested_command"] == "gside db setup --tier sourmash"
        assert method["estimated_size"] == "~3.9GB"
        assert method["error"].startswith("sourmash GTDB database not found")

    def test_mash_missing_db_gets_structured_fields(self, tmp_path):
        payload = _identify_one("q.fna", "mash_refseq", db_dir=str(tmp_path))
        method = payload["methods"]["mash_refseq"]
        assert method["missing_tier"] == "mash"
        assert method["suggested_command"] == "gside db setup --tier mash"
        assert method["estimated_size"] == "~179MB"

    def test_plain_error_has_no_structured_fields(self, monkeypatch):
        import gside.analysis.multigene_identifier as mg

        def boom(c, threads=4):
            raise OSError("blast exploded")

        monkeypatch.setattr(mg, "identify_multigene", boom)
        payload = cli._identify_one("x.fna", "marker", None)
        method = payload["methods"]["marker"]
        assert method["error"] == "blast exploded"
        assert "missing_tier" not in method
        assert "suggested_command" not in method
        assert "estimated_size" not in method

    def test_helper_extracts_known_tiers(self):
        from gside.cli import _structured_missing_db_fields as fields

        assert fields("panel db gone (run: gside db setup --tier panel)") == {
            "suggested_command": "gside db setup --tier panel",
            "missing_tier": "panel",
            "estimated_size": "~130MB sketch / ~1GB genomes",
        }
        assert fields("no setup hint here") == {}
        assert fields("(run: gside db setup --tier )") == {}


class TestAutoSetup:
    def test_flag_parses_and_defaults_off(self, monkeypatch):
        seen = []
        monkeypatch.setattr(cli, "_run_species", lambda *a: seen.append(a) or 0)
        assert CliRunner().invoke(main, ["species", "--auto-setup", "x.fna"]).exit_code == 0
        assert CliRunner().invoke(main, ["species", "x.fna"]).exit_code == 0
        assert seen[0][8] is True
        assert seen[1][8] is False

    def test_missing_tier_triggers_setup_and_retry(self, monkeypatch, capsys):
        calls, setups = [], []

        def fake_identify(contigs, mode, db_dir, threads=4):
            calls.append(mode)
            if len(calls) == 1:
                return _missing_db_payload(contigs)
            return _method_payload(contigs, species="Retried", conf="high")

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: setups.append(tier))
        assert _run_species(["q.fna"], "sourmash", None, auto_setup=True) == 0
        payload = json.loads(capsys.readouterr().out)
        assert setups == ["sourmash"]
        assert calls == ["sourmash", "sourmash"]
        assert payload["methods"]["sourmash"]["result"]["species"] == "Retried"
        assert payload["verdict"]["species"] == "Retried"

    def test_disabled_leaves_error_untouched(self, monkeypatch, capsys):
        setups = []
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _missing_db_payload(c))
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: setups.append(tier))
        assert _run_species(["q.fna"], "sourmash", None) == 0
        payload = json.loads(capsys.readouterr().out)
        assert setups == []
        assert "missing_tier" in payload["methods"]["sourmash"]
        assert payload["verdict"]["species"] == "Unknown"

    def test_retries_only_once(self, monkeypatch, capsys):
        calls, setups = [], []

        def fake_identify(contigs, mode, db_dir, threads=4):
            calls.append(mode)
            return _missing_db_payload(contigs, method="mash_refseq", tier="mash", size="~179MB")

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: setups.append(tier))
        assert _run_species(["q.fna"], "mash_refseq", None, auto_setup=True) == 0
        payload = json.loads(capsys.readouterr().out)
        assert calls == ["mash_refseq", "mash_refseq"]
        assert setups == ["mash"]
        assert payload["methods"]["mash_refseq"]["missing_tier"] == "mash"

    def test_no_missing_tier_no_setup(self, monkeypatch, capsys):
        setups = []
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: setups.append(tier))
        assert _run_species(["q.fna"], "marker", None, auto_setup=True) == 0
        json.loads(capsys.readouterr().out)
        assert setups == []

    def test_setup_failure_keeps_error_entry(self, monkeypatch, capsys):
        calls = []

        def fake_identify(contigs, mode, db_dir, threads=4):
            calls.append(mode)
            return _missing_db_payload(contigs)

        def boom(tier):
            raise RuntimeError("network down")

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", boom)
        assert _run_species(["q.fna"], "sourmash", None, auto_setup=True) == 0
        payload = json.loads(capsys.readouterr().out)
        assert calls == ["sourmash"]
        assert "missing_tier" in payload["methods"]["sourmash"]
        assert payload["verdict"]["species"] == "Unknown"

    def test_setup_progress_goes_to_stderr_not_stdout(self, monkeypatch, capsys):
        import gside.db as dbmod

        calls = []

        def fake_identify(contigs, mode, db_dir, threads=4):
            calls.append(mode)
            if len(calls) == 1:
                return _missing_db_payload(contigs)
            return _method_payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(
            dbmod,
            "db_setup",
            lambda tier, source="": print(f"downloading {tier}...") or {"x": "ok"},
        )
        assert _run_species(["q.fna"], "sourmash", None, auto_setup=True) == 0
        captured = capsys.readouterr()
        assert "downloading sourmash" in captured.err
        assert "downloading sourmash" not in captured.out
        assert json.loads(captured.out)["methods"]["sourmash"]["result"]["species"] == "S"

    def test_all_mode_verdict_recomputed_after_retry(self, monkeypatch, capsys):
        setups = []

        def fake_identify(contigs, mode, db_dir, threads=4):
            if mode == "all":
                payload = _payload(contigs, species="M", conf="high")
                payload["methods"]["sourmash"] = _missing_db_payload(contigs)["methods"]["sourmash"]
                payload["verdict"] = {"species": "M", "confidence": "high", "basis": ["marker"]}
                return payload
            return _method_payload(contigs, species="S", conf="high")

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: setups.append(tier))
        assert _run_species(["q.fna"], "all", None, auto_setup=True) == 0
        payload = json.loads(capsys.readouterr().out)
        assert setups == ["sourmash"]
        assert payload["methods"]["sourmash"]["result"]["species"] == "S"
        assert payload["verdict"]["basis"] == ["sourmash"]
        assert payload["verdict"]["species"] == "S"


class TestGtdbtkMode:
    def test_modes_list_order(self):
        from gside.cli import MODES

        assert MODES == ["marker", "panel", "mash_refseq", "sourmash", "gtdbtk", "all"]

    def test_cli_accepts_gtdbtk_mode(self, monkeypatch):
        seen = {}
        monkeypatch.setattr(cli, "_run_species", lambda *a: seen.update(mode=a[1]) or 0)
        result = CliRunner().invoke(main, ["species", "-m", "gtdbtk", "x.fna"])
        assert result.exit_code == 0
        assert seen["mode"] == "gtdbtk"

    def test_identify_one_dispatches_gtdbtk_with_threads(self, monkeypatch):
        import gside.analysis.gtdbtk_identifier as gmod

        got = {}
        monkeypatch.setattr(
            gmod,
            "identify_by_gtdbtk",
            lambda c, threads=4: got.update(threads=threads) or _FakeResult("G", "high"),
        )
        payload = cli._identify_one("x.fna", "gtdbtk", None, threads=7)
        assert payload["methods"]["gtdbtk"]["species"] == "G"
        assert got["threads"] == 7

    def test_arbitrate_layer3_beats_ani_high(self):
        from gside.cli import _arbitrate

        verdict = _arbitrate(
            {
                "marker": {"species": "M", "confidence": "high"},
                "panel": {"result": {"species": "P", "confidence": "high"}},
                "gtdbtk": {"result": {"species": "G", "confidence": "medium"}},
            }
        )
        assert verdict["basis"] == ["gtdbtk"]
        assert verdict["species"] == "G"

    def test_help_lists_gtdbtk_choice_and_auto_setup(self):
        result = CliRunner().invoke(main, ["scan", "-h"])
        assert result.exit_code == 0
        assert "gtdbtk" in result.output
        assert "--auto-setup" in result.output
        assert "Install a missing database tier" in result.output


class TestCompositeIdentify:
    def test_runs_only_selected_methods(self, monkeypatch):
        import types

        import gside.analysis.multigene_identifier as mg
        import gside.analysis.sourmash_identifier as sm

        def fake_marker(contigs, threads=4):
            return types.SimpleNamespace(
                to_dict=lambda: {
                    "species": "Salmonella",
                    "confidence": "high",
                    "database": {},
                }
            )

        def fake_sourmash(contigs, db_dir=None):
            return types.SimpleNamespace(
                to_dict=lambda: {
                    "species": "Salmonella",
                    "confidence": "medium",
                    "database": {},
                }
            )

        monkeypatch.setattr(mg, "identify_multigene", fake_marker)
        monkeypatch.setattr(sm, "identify_by_sourmash", fake_sourmash)
        payload = _identify_one("g.fna", "marker+sourmash", None)
        assert set(payload["methods"]) == {"marker", "sourmash"}
        assert payload["verdict"]["species"] == "Salmonella"
        assert payload["verdict"]["agreement"] == "2/2"

    def test_deduplicated_modes_run_once(self, monkeypatch):
        import types

        import gside.analysis.multigene_identifier as mg

        calls: list[str] = []

        def fake_marker(contigs, threads=4):
            calls.append("marker")
            return types.SimpleNamespace(
                to_dict=lambda: {"species": "Salmonella", "confidence": "high", "database": {}}
            )

        monkeypatch.setattr(mg, "identify_multigene", fake_marker)
        payload = _identify_one("g.fna", "marker+marker", None)
        assert calls == ["marker"]
        assert list(payload["methods"]) == ["marker"]


class TestCompositeUnknownMode:
    def test_identify_one_rejects_unknown_mode(self):
        with pytest.raises(ValueError, match="unsupported isolate mode"):
            _identify_one("g.fna", "bogus", None)


class TestValidateContract:
    def test_help_lists_output_and_threads(self):
        result = CliRunner().invoke(main, ["validate", "-h"])
        assert result.exit_code == 0
        assert "-o, --output-dir" in result.output
        assert "-t, --threads" in result.output
        assert "--mode" not in result.output

    def test_output_dir_is_required(self, tmp_path):
        source = tmp_path / "g.fna"
        source.write_text(">c\nACGT\n")
        result = CliRunner().invoke(main, ["validate", str(source)])
        assert result.exit_code == 2
        assert "-o" in result.output and "Missing option" in result.output

    def test_mode_flag_is_gone(self, tmp_path):
        source = tmp_path / "g.fna"
        source.write_text(">c\nACGT\n")
        result = CliRunner().invoke(
            main, ["validate", str(source), "-m", "simple", "-o", str(tmp_path / "t")]
        )
        assert result.exit_code == 2
        assert "No such option" in result.output

    def test_output_dir_cannot_be_a_file(self, tmp_path):
        source = tmp_path / "g.fna"
        source.write_text(">c\nACGT\n")
        target = tmp_path / "notadir"
        target.write_text("x")
        result = CliRunner().invoke(main, ["validate", str(source), "-o", str(target)])
        assert result.exit_code == 2
        assert "is a file" in result.output

    def test_partial_status_exits_1(self, monkeypatch, tmp_path):
        import gside.analysis.taxonomic_validator as tv

        monkeypatch.setattr(
            tv,
            "validate_genome",
            lambda c, o, threads=4: tv.ValidateResult(status="partial", errors=["no gtdbtk"]),
        )
        source = tmp_path / "x.fna"
        source.write_text(">c\nACGT\n")
        result = CliRunner().invoke(main, ["validate", str(source), "-o", str(tmp_path / "tax")])
        assert result.exit_code == 1
        assert '"status": "partial"' in result.output
