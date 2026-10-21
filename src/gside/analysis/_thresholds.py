"""Shared ANI confidence thresholds (single source for panel / sylph).

Panel skani and sylph both judge species-level ANI the same way: >= 95 is a
high-confidence same-species call, 93-95 is medium, below 93 is low.
"""

from __future__ import annotations

ANI_HIGH = 95.0
ANI_MEDIUM = 93.0


def ani_confidence(ani: float) -> str:
    if ani >= ANI_HIGH:
        return "high"
    if ani >= ANI_MEDIUM:
        return "medium"
    return "low"
