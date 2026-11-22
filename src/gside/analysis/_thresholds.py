"""Shared ANI confidence thresholds (single source for panel engines).

Panel engines (skani/sourmash) judge species-level ANI the same way: >= 95 is a
high-confidence same-species call, 93-95 is medium, below 93 is low.
"""

from __future__ import annotations

from enum import StrEnum

ANI_HIGH = 95.0
ANI_MEDIUM = 93.0


class Confidence(StrEnum):
    """Canonical species-identification confidence tiers."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


def ani_confidence(ani: float) -> Confidence:
    if ani >= ANI_HIGH:
        return Confidence.HIGH
    if ani >= ANI_MEDIUM:
        return Confidence.MEDIUM
    return Confidence.LOW
