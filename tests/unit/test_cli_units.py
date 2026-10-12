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
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d: _payload(c))
        assert _run_species(["a.fna"], "marker", None) == 0
        assert json.loads(capsys.readouterr().out)["contigs"] == "a.fna"

    def test_json_multi_is_array(self, monkeypatch, capsys):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d: _payload(c))
        assert _run_species(["a.fna", "b.fna"], "marker", None) == 0
        docs = json.loads(capsys.readouterr().out)
        assert [d["contigs"] for d in docs] == ["a.fna", "b.fna"]

    def test_tsv_and_md(self, monkeypatch, capsys):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d: _payload(c))
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
            lambda contigs: [
                {"gene": "inva", "identity": 80.0, "coverage": 50.0},
            ],
        )
        r = mg.identify_multigene("dummy.fna").to_dict()
        assert r["species"] == "Unknown"
        assert any("Near-threshold" in n for n in r["notes"])

    def test_no_hits_no_notes(self, monkeypatch):
        import gside.analysis.multigene_identifier as mg

        monkeypatch.setattr(mg, "_blast_contigs", lambda contigs: [])
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
        assert "markers" in result.output

    def test_list(self, tmp_path, monkeypatch, capsys):
        from click.testing import CliRunner

        self._db(tmp_path, monkeypatch)
        result = CliRunner().invoke(main, ["db", "list"])
        assert result.exit_code == 0
        assert "sourmash" in result.output

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
        import gside.analysis.multigene_identifier as mg
        import gside.analysis.sourmash_identifier as smod
        import gside.cli as climod

        monkeypatch.setattr(mg, "identify_multigene", lambda c: _FakeResult("M", "high"))
        monkeypatch.setattr(
            ani, "identify_by_ani", lambda c, mode, db_dir: _FakeResult("A", "high")
        )
        monkeypatch.setattr(
            smod, "identify_by_sourmash", lambda c, db_dir: _FakeResult("S", "high")
        )
        payload = climod._identify_one("x.fna", "all", None)
        assert set(payload["methods"]) == {"marker", "panel", "mash_refseq", "sourmash"}
        assert payload["verdict"]["basis"] == ["panel"]

    def test_branch_error_entry(self, monkeypatch):
        import gside.analysis.multigene_identifier as mg
        import gside.cli as climod

        def boom(c):
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

        class _R:
            def to_dict(self):
                return {"ok": True}

        monkeypatch.setattr(
            "gside.analysis.taxonomic_validator.validate_genome",
            lambda c, mode, output_dir: _R(),
        )
        result = CliRunner().invoke(main, ["validate", "x.fna", "--mode", "simple"])
        assert result.exit_code == 0 and '"ok": true' in result.output.lower()


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
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d: _payload(c))
        out = tmp_path / "out.json"
        result = CliRunner().invoke(main, ["species", "--output", str(out), "a.fna"])
        assert result.exit_code == 0
        assert json.loads(out.read_text())["contigs"] == "a.fna"
        assert result.stdout == ""
        assert f"Results written to {out}" in result.stderr

    def test_output_short_flag_multi_file_array(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d: _payload(c))
        out = tmp_path / "o.json"
        result = CliRunner().invoke(main, ["species", "-o", str(out), "a.fna", "b.fna"])
        assert result.exit_code == 0
        assert [d["contigs"] for d in json.loads(out.read_text())] == ["a.fna", "b.fna"]
        assert "Results written to" in result.stderr

    def test_output_tsv_payload(self, monkeypatch, tmp_path):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d: _payload(c))
        out = tmp_path / "o.tsv"
        result = CliRunner().invoke(main, ["species", "-f", "tsv", "-o", str(out), "a.fna"])
        assert result.exit_code == 0
        assert out.read_text().startswith("contigs\t")
        assert result.stdout == ""


class TestDbOutputFormats:
    _STATUS = {
        "markers": {"ready": True, "path": "/p/D1_marker", "tier": "mini"},
        "D2_ani": {"ready": False, "path": "/p/D2_ani", "tier": "panel"},
    }

    _TIERS = {
        "mini": {
            "description": "Marker rules + sequences (bundled)",
            "size": "~1MB",
            "items": ["markers"],
        },
    }

    def test_status_tsv_machine_readable(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "status", "--format", "tsv"])
        assert result.exit_code == 0
        lines = result.output.splitlines()
        assert lines[0] == "COMPONENT\tREADY\tTIER\tPATH"
        assert lines[1] == "markers\tready\tmini\t/p/D1_marker"
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

    def test_status_table_default_keeps_icons_and_hint(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: self._STATUS)
        result = CliRunner().invoke(main, ["db", "status"])
        assert result.exit_code == 0
        assert "✅" in result.output and "❌" in result.output
        assert "Run 'gside db setup --tier <tier>' to provision" in result.output

    def test_list_tsv(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "_TIERS", self._TIERS)
        result = CliRunner().invoke(main, ["db", "list", "--format", "tsv"])
        assert result.exit_code == 0
        assert result.output.splitlines() == [
            "NAME\tSIZE\tDESCRIPTION",
            "mini\t~1MB\tMarker rules + sequences (bundled)",
        ]

    def test_list_json(self, monkeypatch):
        import gside.db as dbmod

        monkeypatch.setattr(dbmod, "_TIERS", self._TIERS)
        result = CliRunner().invoke(main, ["db", "list", "-f", "json"])
        assert result.exit_code == 0
        assert json.loads(result.output) == self._TIERS


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
        assert v == {"species": "Unknown", "confidence": "low", "basis": ["a", "b"]}

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
