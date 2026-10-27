"""Shared fixtures for hermetic unit tests.

Most unit tests pin the legacy subprocess implementation, which needs no gseek
import. Tests that exercise gseek or shadow set ``GSIDE_ENGINE`` explicitly;
the default-mode test deletes the variable to assert the shipped default.
"""

from __future__ import annotations

import os
from collections.abc import Iterator

import pytest


@pytest.fixture(autouse=True)
def _legacy_engine_for_hermetic_tests() -> Iterator[None]:
    previous = os.environ.get("GSIDE_ENGINE")
    os.environ["GSIDE_ENGINE"] = "legacy"
    yield
    if previous is None:
        os.environ.pop("GSIDE_ENGINE", None)
    else:
        os.environ["GSIDE_ENGINE"] = previous
