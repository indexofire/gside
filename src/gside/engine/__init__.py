"""Alignment backends for gside's species-identification hot paths.

Live consumers: :class:`~gside.engine.backends.blast.BlastBackend` (marker),
:class:`~gside.engine.backends.skani.SkaniBackend` (panel), and
:class:`~gside.engine.backends.kmer.MashBackend` / ``SourmashBackend``
(mash_refseq / sourmash). Binary discovery and checked subprocess execution
live in :mod:`gside.engine._env`; :class:`gside.engine.hits.Hit` is the shared
BLAST row representation.
"""

from __future__ import annotations

from .hits import Hit

__all__ = ["Hit"]
