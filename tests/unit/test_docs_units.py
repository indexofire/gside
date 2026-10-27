"""Docs-consistency tests: user docs must match the CLI surface."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import yaml

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

    def test_validate_documents_checkm2_and_gtdbtk(self):
        en, zh = _read("usage/cli.md"), _read("usage/cli.zh.md")
        for token in ("CheckM2", "GTDB-Tk", "validation.json", "required"):
            assert token in en, token
        assert "CheckM2" in zh and "GTDB-Tk" in zh and "validation.json" in zh

    def test_species_auto_setup_documented(self):
        en, zh = _read("usage/cli.md"), _read("usage/cli.zh.md")
        assert "--auto-setup" in en
        assert "--auto-setup" in zh


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
        assert "data/db/D1_marker/marker_rules.yaml" in en


class TestBilingualDocs:
    def _i18n(self) -> dict:
        config = yaml.safe_load((_REPO / "mkdocs.yml").read_text(encoding="utf-8"))
        for plugin in config["plugins"]:
            if isinstance(plugin, dict) and "i18n" in plugin:
                return plugin["i18n"]
        raise AssertionError("mkdocs.yml has no i18n plugin")

    def test_language_selector_is_generated_by_static_i18n(self):
        i18n = self._i18n()
        assert i18n["reconfigure_material"] is True
        assert i18n["reconfigure_search"] is True

    def test_chinese_navigation_titles_are_translated(self):
        languages = {lang["locale"]: lang for lang in self._i18n()["languages"]}
        translations = languages["zh"]["nav_translations"]
        assert translations["Home"] == "首页"
        assert translations["Usage"] == "使用"
        assert translations["Development"] == "开发"

    def test_every_english_page_has_a_chinese_pair(self):
        missing = [
            str(page.relative_to(_DOCS))
            for page in _DOCS.rglob("*.md")
            if not page.name.endswith(".zh.md")
            and not page.with_name(f"{page.stem}.zh.md").exists()
        ]
        assert not missing


class TestGseekFloor:
    def test_pyproject_and_pixi_floors_agree(self):
        pyproject = tomllib.loads((_REPO / "pyproject.toml").read_text(encoding="utf-8"))
        floor_pyproject = next(
            dep.split(">=", 1)[1]
            for dep in pyproject["project"]["dependencies"]
            if dep.startswith("gseek>=")
        )
        pixi = (_REPO / "pixi.toml").read_text(encoding="utf-8")
        floor_pixi = re.search(r'^gseek = ">=([^"]+)"', pixi, re.MULTILINE).group(1)
        assert floor_pyproject == floor_pixi
        assert tuple(map(int, floor_pyproject.split("."))) >= (0, 10, 2)


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
