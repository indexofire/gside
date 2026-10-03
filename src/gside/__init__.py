"""gside — Genome Species IDentification Engine."""

from __future__ import annotations

try:
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version("gside")
    except PackageNotFoundError:
        __version__ = "0.1.0.dev"
except Exception:  # pragma: no cover — extremely defensive fallback
    __version__ = "0.1.0.dev"

__all__ = ["__version__"]
