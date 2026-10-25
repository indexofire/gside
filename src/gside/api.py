"""Stable, import-friendly access to gside's identification engines.

The CLI adds database provisioning, layered arbitration, and output formatting;
these functions expose the individual engines for notebooks, pipelines, and
tests. They require the same external tools and reference databases as the
corresponding CLI mode.

The API is intentionally small and may gain optional arguments before 1.0;
existing positional arguments and return types are kept stable.
"""

from __future__ import annotations

from gside.analysis.ani_identifier import AniIdResult, identify_by_ani
from gside.analysis.gtdbtk_identifier import GtdbkResult, identify_by_gtdbtk
from gside.analysis.kraken2_identifier import Kraken2Result, identify_by_kraken2
from gside.analysis.multigene_identifier import MultiGeneResult, identify_multigene
from gside.analysis.read_screening import screen_reads
from gside.analysis.sourmash_identifier import SourmashIdResult, identify_by_sourmash
from gside.analysis.sylph_identifier import SylphResult, identify_by_sylph
from gside.analysis.taxonomic_validator import ValidateResult, validate_genome
from gside.marker_db import MarkerDatabase

__all__ = [
    "AniIdResult",
    "GtdbkResult",
    "Kraken2Result",
    "MarkerDatabase",
    "MultiGeneResult",
    "SourmashIdResult",
    "SylphResult",
    "ValidateResult",
    "identify_by_ani",
    "identify_by_gtdbtk",
    "identify_by_kraken2",
    "identify_by_sourmash",
    "identify_by_sylph",
    "identify_multigene",
    "screen_reads",
    "validate_genome",
]
