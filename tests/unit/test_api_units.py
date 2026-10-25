"""Public Python API surface (no external tools)."""

from __future__ import annotations

import gside.api as api


def test_api_exports_are_importable() -> None:
    for name in api.__all__:
        assert hasattr(api, name), name


def test_api_all_is_unique() -> None:
    assert len(api.__all__) == len(set(api.__all__))


def test_api_covers_each_engine_family() -> None:
    assert "identify_multigene" in api.__all__
    assert "identify_by_ani" in api.__all__
    assert "identify_by_sourmash" in api.__all__
    assert "identify_by_gtdbtk" in api.__all__
    assert "identify_by_kraken2" in api.__all__
    assert "identify_by_sylph" in api.__all__
    assert "screen_reads" in api.__all__
    assert "validate_genome" in api.__all__
