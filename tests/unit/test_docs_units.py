"""Docs-consistency tests: user docs must match the CLI surface."""

from __future__ import annotations

from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
_DOCS = _REPO / "docs"


def _read(*parts: str) -> str:
    return (_DOCS.joinpath(*parts)).read_text(encoding="utf-8")


class TestCliModesDocumented:
    def test_species_modes(self):
        from gside.cli import MODES

        en, zh = _read("usage/cli.md"), _read("usage/cli.zh.md")
        for mode in MODES:
            assert f"`{mode}`" in en, mode
            assert f"`{mode}`" in zh, mode

    def test_db_tiers(self):
        from gside.cli import TIERS

        en, zh = _read("usage/cli.md"), _read("usage/cli.zh.md")
        for tier in TIERS:
            assert tier in en, tier
            assert tier in zh, tier

    def test_validate_modes(self):
        en, zh = _read("usage/cli.md"), _read("usage/cli.zh.md")
        for mode in ("simple", "standard"):
            assert mode in en and mode in zh


class TestDatabaseDocs:
    def test_tier_dirs(self):
        from gside import db as dbmod

        en, zh = _read("installation/databases.md"), _read("installation/databases.zh.md")
        for tier, info in dbmod._TIERS.items():
            for item in info["items"]:
                assert item in en, (tier, item)
                assert item in zh, (tier, item)

    def test_marker_paths(self):
        en = _read("reference/marker-rules.md")
        assert "data/db/L1_marker/marker_rules.yaml" in en


class TestCoverageFloorConsistency:
    def _floor_sources(self):
        import tomllib

        repo = _REPO
        with open(repo / "pyproject.toml", "rb") as fh:
            pyproject = tomllib.load(fh)
        floor_pyproject = pyproject["tool"]["coverage"]["report"]["fail_under"]
        pixi = (repo / "pixi.toml").read_text()
        import re

        floor_pixi = int(re.search(r"cov-fail-under=(\d+)", pixi).group(1))
        ci = (repo / ".github" / "workflows" / "test.yml").read_text()
        floor_ci = int(re.search(r"cov-fail-under=(\d+)", ci).group(1))
        return floor_pyproject, floor_pixi, floor_ci

    def test_floors_agree(self):
        assert len(set(self._floor_sources())) == 1
