"""Pure unit tests for config env helpers and package version fallback."""

from __future__ import annotations

from gside.config import _env_path


class TestEnvPath:
    def test_first_hit_wins(self, monkeypatch, tmp_path):
        monkeypatch.setenv("GSIDE_TEST_A", "")
        monkeypatch.setenv("GSIDE_TEST_B", str(tmp_path))
        assert _env_path("GSIDE_TEST_A", "GSIDE_TEST_B") == tmp_path

    def test_all_missing_is_none(self, monkeypatch):
        monkeypatch.delenv("GSIDE_TEST_MISSING", raising=False)
        assert _env_path("GSIDE_TEST_MISSING") is None

    def test_checkm2_gtdb_defaults_unset(self, monkeypatch):
        import gside.config as config_mod

        monkeypatch.delenv("CHECKM2DB", raising=False)
        monkeypatch.delenv("GTDBTK_DATA_PATH", raising=False)
        monkeypatch.delenv("GTDBDB", raising=False)
        import importlib

        reloaded = importlib.reload(config_mod)
        assert reloaded.CHECKM2_DB is None and reloaded.GTDB_DB is None

    def test_gtdbtk_data_path_preferred(self, monkeypatch, tmp_path):
        import gside.config as config_mod

        monkeypatch.setenv("GTDBTK_DATA_PATH", str(tmp_path / "a"))
        monkeypatch.setenv("GTDBDB", str(tmp_path / "b"))
        assert config_mod._env_path("GTDBTK_DATA_PATH", "GTDBDB") == tmp_path / "a"


class TestVersionFallback:
    def test_missing_dist_info(self, monkeypatch):
        import importlib
        import importlib.metadata

        def boom(name):
            raise importlib.metadata.PackageNotFoundError

        monkeypatch.setattr(importlib.metadata, "version", boom)
        import gside

        assert importlib.reload(gside).__version__ == "0.1.0.dev"
        monkeypatch.undo()
        importlib.reload(gside)
