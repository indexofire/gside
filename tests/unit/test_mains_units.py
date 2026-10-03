"""Entry-point guard tests (runpy, --help only, no binaries run)."""

from __future__ import annotations

import runpy

import pytest


def _run_main_help(modname, monkeypatch):
    monkeypatch.setattr("sys.argv", ["prog", "--help"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_module(modname, run_name="__main__", alter_sys=True)
    assert exc.value.code == 0


class TestMainGuards:
    def test_package_main(self, monkeypatch):
        _run_main_help("gside", monkeypatch)

    def test_ani_main(self, monkeypatch):
        _run_main_help("gside.analysis.ani_identifier", monkeypatch)

    def test_multigene_main(self, monkeypatch):
        _run_main_help("gside.analysis.multigene_identifier", monkeypatch)

    def test_cli_module_main(self, monkeypatch):
        _run_main_help("gside.cli", monkeypatch)
