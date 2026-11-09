"""Named per-mode database resolution (pure, no binaries)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gside.db_registry import (
    DbChoice,
    enumerate_databases,
    layout_valid,
    parse_db_value,
    resolve_db,
    resolve_selections,
)


def _marker_db(root: Path, name: str = "food") -> Path:
    path = root / "marker" / name
    path.mkdir(parents=True)
    (path / "markers.fasta").write_text(">markers~~~g1~~~x marker\nACGTACGT\n")
    (path / "marker_rules.yaml").write_text("rules:\n- species: X\n  genes: [g1]\n")
    return path


def _panel_db(root: Path, name: str = "lab2026") -> Path:
    path = root / "panel" / name
    (path / "panel.sketch" / "sketches.db").parent.mkdir(parents=True)
    (path / "panel.sketch" / "sketches.db").write_bytes(b"sketch")
    (path / "metadata.tsv").write_text("fna\tspecies\n")
    return path


class TestParseDbValue:
    def test_mode_qualified(self):
        assert parse_db_value("panel:lab2026") == ("panel", "lab2026")

    def test_bare_name(self):
        assert parse_db_value("food") == (None, "food")

    def test_absolute_path_is_bare(self):
        mode, value = parse_db_value("/opt/dbs/food")
        assert mode is None
        assert value == "/opt/dbs/food"

    def test_empty_name_rejected(self):
        with pytest.raises(ValueError, match="empty database name"):
            parse_db_value("panel:")

    def test_empty_value_rejected(self):
        with pytest.raises(ValueError, match="empty database name"):
            parse_db_value("")


class TestResolveDb:
    def test_existing_path_is_used_directly(self, tmp_path):
        db = _marker_db(tmp_path)
        choice = resolve_db("marker", str(db))
        assert choice.mode == "marker"
        assert choice.name == "food"
        assert choice.path == db

    def test_named_under_family_dir(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        db = _marker_db(tmp_path)
        choice = resolve_db("marker", "food")
        assert choice.path == db
        assert choice.name == "food"

    def test_default_selects_the_bundled_tier(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        choice = resolve_db("panel", "default")
        assert choice.name == "default"
        assert choice.path == tmp_path / "D2_ani"

    def test_unknown_name_lists_available(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        _marker_db(tmp_path)
        with pytest.raises(ValueError, match="no marker database named 'nope'"):
            resolve_db("marker", "nope")

    def test_invalid_layout_rejected(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        (tmp_path / "marker" / "broken").mkdir(parents=True)
        with pytest.raises(ValueError, match="not a valid marker database"):
            resolve_db("marker", "broken")

    def test_mode_without_family_rejected(self):
        with pytest.raises(ValueError, match="gtdbtk does not support named databases"):
            resolve_db("gtdbtk", "whatever")

    def test_root_override(self, tmp_path):
        db = _marker_db(tmp_path / "custom-root")
        choice = resolve_db("marker", "food", root=tmp_path / "custom-root")
        assert choice.path == db


class TestResolveSelections:
    def test_qualified_per_mode(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        marker = _marker_db(tmp_path)
        panel = _panel_db(tmp_path)
        selections = resolve_selections(["marker:food", "panel:lab2026"], ["marker", "panel"])
        assert [c.path for c in selections["marker"]] == [marker]
        assert [c.path for c in selections["panel"]] == [panel]

    def test_bare_name_needs_single_db_mode(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        marker = _marker_db(tmp_path)
        assert [c.path for c in resolve_selections(["food"], ["marker"])["marker"]] == [marker]

    def test_bare_name_with_multiple_modes_rejected(self):
        with pytest.raises(ValueError, match="qualify the database"):
            resolve_selections(["food"], ["marker", "panel"])

    def test_unknown_mode_prefix_rejected(self):
        with pytest.raises(ValueError, match="not part of this invocation"):
            resolve_selections(["kraken2:std"], ["marker"])

    def test_no_family_mode_rejected(self):
        with pytest.raises(ValueError, match="no mode in this invocation supports"):
            resolve_selections(["std"], ["gtdbtk"])

    def test_duplicates_deduplicated_in_order(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        marker = _marker_db(tmp_path)
        selections = resolve_selections(["food", "marker:food"], ["marker"])
        assert len(selections["marker"]) == 1
        assert selections["marker"][0].path == marker

    def test_mixed_path_and_name(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        named = _marker_db(tmp_path)
        other = _marker_db(tmp_path, name="extra")
        selections = resolve_selections([str(other), "food"], ["marker"])
        assert [c.path for c in selections["marker"]] == [other, named]


class TestEnumerate:
    def test_lists_default_and_named(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        _marker_db(tmp_path)
        bundled = tmp_path / "D1_marker"
        bundled.mkdir()
        (bundled / "markers.fasta").write_text(">markers~~~g1~~~x marker\nACGT\n")
        (bundled / "marker_rules.yaml").write_text("rules:\n- species: X\n  genes: [g1]\n")
        databases = enumerate_databases()
        names = [(c.name, c.ready) for c in databases["marker"]]
        assert ("default", True) in names
        assert ("food", True) in names

    def test_default_not_provisioned_is_not_ready(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        databases = enumerate_databases()
        default = next(c for c in databases["panel"] if c.name == "default")
        assert default.ready is False

    def test_invalid_entries_skipped(self, tmp_path, monkeypatch):
        from gside import config

        monkeypatch.setattr(config, "SPECIES_DB_DIR", tmp_path)
        (tmp_path / "panel" / "junk").mkdir(parents=True)
        databases = enumerate_databases()
        assert all(c.name != "junk" for c in databases["panel"])


def test_db_choice_shape(tmp_path):
    choice = DbChoice(mode="marker", name="food", path=tmp_path, ready=True)
    assert choice.mode == "marker"


class TestLayoutValidation:
    def test_mash_layout(self, tmp_path):
        (tmp_path / "mash.msh").write_bytes(b"msh")
        assert layout_valid("mash_refseq", tmp_path) is True
        assert layout_valid("panel", tmp_path) is False

    def test_sourmash_layout_needs_both_files(self, tmp_path):
        assert layout_valid("sourmash", tmp_path) is False
        (tmp_path / "gtdb-reps-k31.zip").write_bytes(b"zip")
        assert layout_valid("sourmash", tmp_path) is False
        (tmp_path / "lineages.csv").write_text("a,b\n")
        assert layout_valid("sourmash", tmp_path) is True

    def test_panel_accepts_sketch_file_variant(self, tmp_path):
        (tmp_path / "panel.sketch").write_bytes(b"sketch-file")
        assert layout_valid("panel", tmp_path) is True

    def test_missing_directory_is_invalid(self, tmp_path):
        assert layout_valid("marker", tmp_path / "nope") is False


class TestParseEdges:
    def test_leading_colon_rejected(self):
        with pytest.raises(ValueError, match="empty database name"):
            parse_db_value(":panel")

    def test_spaces_are_stripped(self):
        assert parse_db_value(" panel : lab2026 ") == ("panel", "lab2026")


def test_layout_requirements_rejects_unknown_mode():
    from gside.db_registry import _layout_requirements

    with pytest.raises(ValueError, match="does not support named databases"):
        _layout_requirements("gtdbtk")


def test_resolve_db_path_with_wrong_layout(tmp_path):
    (tmp_path / "not-a-panel").mkdir()
    with pytest.raises(ValueError, match="is not a valid panel database"):
        resolve_db("panel", str(tmp_path / "not-a-panel"))
