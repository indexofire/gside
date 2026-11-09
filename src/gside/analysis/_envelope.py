"""Shared JSON envelope for per-method result dataclasses.

Every identifier returns ``{analysis_type, method, database, result}``; the
base class below is the single home for that shape. Subclasses only override
the class-level ``analysis_type`` / ``method`` defaults.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar

SPECIES = "species_identification"
METAGENOMIC = "metagenomic_profiling"


@dataclass
class MethodEnvelope:
    analysis_type: ClassVar[str] = SPECIES
    method: str = ""
    database: dict[str, str] = field(default_factory=dict)
    result: dict[str, Any] = field(default_factory=dict)
    databases: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "analysis_type": self.analysis_type,
            "method": self.method,
            "database": self.database,
            "result": self.result,
        }
        if self.databases:
            payload["databases"] = self.databases
        return payload
