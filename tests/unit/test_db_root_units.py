"""Database-root resolution, root reporting, and bundled-D1 seeding (pure).

Pins the gspec-style datadir chain (env > repo checkout > XDG data home),
`db status` root reporting (list stays silent), and the wheel-bundled
D1_marker seed materialized into a fresh user root.
"""

from __future__ import annotations

import json
from pathlib import Path

from gside import db
from gside.config import (
    DB_ROOT_SOURCE,
    SPECIES_DB_DIR,
    DbRootSource,
    db_root_source,
    resolve_db_root,
)

_NO_REPO = Path("/nonexistent-repo-data")


def _fake_bundle(tmp_path: Path) -> Path:
    """A wheel-like BUNDLED_DB_ROOT carrying D1_marker/default + cpin."""
    bundled = tmp_path / "wheel" / "data" / "db"
    for name, gene in (("default", "g1"), ("cpin", "c1")):
        entry = bundled / "D1_marker" / name
        entry.mkdir(parents=True)
        (entry / "markers.fasta").write_text(
            f">gside|db=test|gene={gene}|acc=t|func=primary\nACGT\n"
        )
        (entry / "marker_rules.yaml").write_text(
            f"rules:\n- species: {name.title()}\n  genes: [{gene}]\n"
        )
    return bundled


class TestResolveDbRoot:
    def test_env_db_dir_wins(self, tmp_path: Path) -> None:
        path, source = resolve_db_root(
            {"GSIDE_DB_DIR": str(tmp_path / "custom-db")}, repo_data_dir=tmp_path
        )
        assert path == tmp_path / "custom-db"
        assert source == DbRootSource.ENV

    def test_env_data_dir_appends_db(self, tmp_path: Path) -> None:
        path, source = resolve_db_root(
            {"GSIDE_DATA_DIR": str(tmp_path / "data")}, repo_data_dir=tmp_path
        )
        assert path == tmp_path / "data" / "db"
        assert source == DbRootSource.ENV

    def test_env_db_dir_beats_data_dir(self, tmp_path: Path) -> None:
        path, _ = resolve_db_root(
            {"GSIDE_DB_DIR": str(tmp_path / "a"), "GSIDE_DATA_DIR": str(tmp_path / "b")},
            repo_data_dir=_NO_REPO,
        )
        assert path == tmp_path / "a"

    def test_blank_env_values_are_ignored(self, tmp_path: Path) -> None:
        path, source = resolve_db_root(
            {"GSIDE_DB_DIR": "  ", "GSIDE_DATA_DIR": ""}, repo_data_dir=tmp_path
        )
        assert source == DbRootSource.REPO
        assert path == tmp_path / "db"

    def test_repo_checkout_wins_over_xdg(self, tmp_path: Path) -> None:
        path, source = resolve_db_root(
            {"XDG_DATA_HOME": str(tmp_path / "xdg")}, repo_data_dir=tmp_path
        )
        assert path == tmp_path / "db"
        assert source == DbRootSource.REPO

    def test_xdg_data_home_absolute(self, tmp_path: Path) -> None:
        path, source = resolve_db_root(
            {"XDG_DATA_HOME": str(tmp_path / "xdg")}, repo_data_dir=_NO_REPO
        )
        assert path == tmp_path / "xdg" / "gside" / "db"
        assert source == DbRootSource.XDG

    def test_xdg_data_home_relative_is_ignored(self) -> None:
        path, source = resolve_db_root({"XDG_DATA_HOME": "relative/path"}, repo_data_dir=_NO_REPO)
        assert path == Path.home() / ".local" / "share" / "gside" / "db"
        assert source == DbRootSource.XDG

    def test_import_time_constants_match_the_resolver(self) -> None:
        assert resolve_db_root()[0] == SPECIES_DB_DIR
        assert resolve_db_root()[1] == DB_ROOT_SOURCE


class TestDbRootSourceClassification:
    def test_env_candidate(self, tmp_path: Path) -> None:
        env = {"GSIDE_DB_DIR": str(tmp_path / "custom-db")}
        assert db_root_source(tmp_path / "custom-db", env, repo_data_dir=_NO_REPO) == "env"

    def test_repo_candidate(self, tmp_path: Path) -> None:
        assert db_root_source(tmp_path / "db", {}, repo_data_dir=tmp_path) == "repo"

    def test_xdg_candidate(self, tmp_path: Path) -> None:
        env = {"XDG_DATA_HOME": str(tmp_path / "xdg")}
        assert (
            db_root_source(tmp_path / "xdg" / "gside" / "db", env, repo_data_dir=_NO_REPO) == "xdg"
        )

    def test_unknown_path_is_custom(self, tmp_path: Path) -> None:
        assert db_root_source(tmp_path, {}, repo_data_dir=_NO_REPO) == "custom"


class TestStatusReportsRoot:
    def test_status_json_carries_root_and_source(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "SPECIES_DB_DIR", tmp_path)
        payload = json.loads(db.format_status("json"))
        assert set(payload) == {"modes", "tiers", "root", "root_source"}
        assert payload["root"] == str(tmp_path)
        assert payload["root_source"] == "custom"
        assert "modes" in payload and "tiers" in payload

    def test_status_table_prepends_root_line(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "SPECIES_DB_DIR", tmp_path)
        text = db.format_status("table")
        first = text.splitlines()[0]
        assert first.startswith("Database root:")
        assert str(tmp_path) in first
        assert "gside databases" in text

    def test_status_tsv_and_csv_unchanged(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "SPECIES_DB_DIR", tmp_path)
        assert "Database root" not in db.format_status("tsv")
        assert "Database root" not in db.format_status("csv")

    def test_list_view_has_no_root_keys(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "SPECIES_DB_DIR", tmp_path)
        for fmt in ("json", "table", "tsv", "csv"):
            assert "root_source" not in db.format_databases(fmt)
            assert "Database root" not in db.format_databases(fmt)


class TestMarkerSeeding:
    def test_seed_copies_bundled_d1_into_root(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        root = tmp_path / "xdg" / "gside" / "db"
        seeded = db.ensure_marker_seed(root)
        assert sorted(seeded) == ["cpin", "default"]
        assert (root / "D1_marker" / "default" / "marker_rules.yaml").is_file()
        assert (root / "D1_marker" / "cpin" / "markers.fasta").is_file()

    def test_seed_skips_stray_files_in_bundled_d1(self, tmp_path: Path, monkeypatch) -> None:
        bundled = _fake_bundle(tmp_path)
        (bundled / "D1_marker" / "README.txt").write_text("stray file\n")
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", bundled)
        root = tmp_path / "root"
        seeded = db.ensure_marker_seed(root)
        assert sorted(seeded) == ["cpin", "default"]
        assert not (root / "D1_marker" / "README.txt").exists()

    def test_seed_creates_a_private_root(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        root = tmp_path / "xdg" / "gside" / "db"
        db.ensure_marker_seed(root)
        assert root.is_dir()
        assert (root.stat().st_mode & 0o700) == 0o700

    def test_seed_is_idempotent(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        root = tmp_path / "root"
        assert db.ensure_marker_seed(root)
        assert db.ensure_marker_seed(root) == []

    def test_seed_without_bundle_creates_nothing(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", tmp_path / "absent")
        root = tmp_path / "root"
        assert db.ensure_marker_seed(root) == []
        assert not root.exists()

    def test_seed_swallows_copy_errors(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))

        def boom(*args, **kwargs):
            raise OSError("read-only filesystem")

        monkeypatch.setattr(db.shutil, "copytree", boom)
        assert db.ensure_marker_seed(tmp_path / "root") == []

    def test_component_ready_falls_back_to_bundled(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "SPECIES_DB_DIR", tmp_path / "empty-root")
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        assert db._component_ready("D1_marker") is True
        assert db._component_ready("D2_ani") is False

    def test_db_setup_mini_seeds_a_fresh_root(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        root = tmp_path / "root"
        monkeypatch.setattr(db, "SPECIES_DB_DIR", root)
        results = db.db_setup("mini")
        assert results["D1_marker"].startswith("materialized from bundled")
        assert (root / "D1_marker" / "default" / "marker_rules.yaml").is_file()

    def test_format_databases_seeds_marker(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        root = tmp_path / "root"
        monkeypatch.setattr(db, "SPECIES_DB_DIR", root)
        db.format_databases("tsv")
        assert (root / "D1_marker" / "default" / "marker_rules.yaml").is_file()


class TestDefaultMarkerDbSeeding:
    def test_seeds_when_root_is_missing(self, tmp_path: Path, monkeypatch) -> None:
        from gside import config
        from gside.marker_db import default_marker_db

        root = tmp_path / "root"
        monkeypatch.setattr(config, "SPECIES_DB_DIR", root)
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        database = default_marker_db()
        assert database.root == root / "D1_marker" / "default"
        assert database.name == "species_markers_v2"
        assert (database.root / "markers.fasta").is_file()


class TestIdentifyMultigeneSeeding:
    def test_default_path_uses_seeded_marker_db(self, tmp_path: Path, monkeypatch) -> None:
        import gside.analysis.multigene_identifier as mg

        monkeypatch.delenv("GSIDE_MARKER_DB", raising=False)
        root = tmp_path / "root"
        monkeypatch.setattr(
            mg, "_MARKER_RULES", root / "D1_marker" / "default" / "marker_rules.yaml"
        )
        monkeypatch.setattr(mg, "_MARKERS_FASTA", root / "D1_marker" / "default" / "markers.fasta")
        monkeypatch.setattr(
            mg, "_MARKER_METADATA", root / "D1_marker" / "default" / "metadata.yaml"
        )
        monkeypatch.setattr(db, "BUNDLED_DB_ROOT", _fake_bundle(tmp_path))
        monkeypatch.setattr(
            mg,
            "_blast_contigs",
            lambda contigs, threads=4, marker_db=None: [
                {"gene": "g1", "identity": 99.0, "coverage": 100.0}
            ],
        )
        result = mg.identify_multigene("x.fna").to_dict()
        assert result["species"] == "Default"
        assert (root / "D1_marker" / "default" / "marker_rules.yaml").is_file()


class TestEffectiveRootReported:
    def test_status_reports_the_effective_root_over_the_wire(self, tmp_path: Path) -> None:
        # db_root_source defaults to the effective SPECIES_DB_DIR; on this
        # machine (source checkout) that is the repo candidate.
        assert db_root_source(SPECIES_DB_DIR) in ("env", "repo", "xdg")
