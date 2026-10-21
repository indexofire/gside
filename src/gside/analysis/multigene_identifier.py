"""Multi-gene species identifier using the expanded markers BLAST database.

Replaces the single-gene species_identifier for routine use. Reads
marker_rules.yaml for gene combination logic and scoring.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from gside.config import SPECIES_DB_DIR

_MARKERS_FASTA = SPECIES_DB_DIR / "D1_marker" / "markers.fasta"
_MARKER_RULES = SPECIES_DB_DIR / "D1_marker" / "marker_rules.yaml"

_MIN_IDENTITY = 85.0
_MIN_COVERAGE = 60.0
_HIGH_CONF = 90.0


@dataclass
class MultiGeneResult:
    species: str = "Unknown"
    confidence: str = "low"
    detected_markers: list[dict[str, Any]] = field(default_factory=list)
    all_hits: list[dict[str, Any]] = field(default_factory=list)
    matched_rule: str = ""
    method: str = "multigene"
    database_version: str = "unknown"
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "confidence": self.confidence,
            "method": self.method,
            "database": {"name": "species_markers_v2", "version": self.database_version},
            "detected_markers": self.detected_markers,
            "matched_rule": self.matched_rule,
            "notes": self.notes,
            "interpretation": self._interpret(),
        }

    def _interpret(self) -> str:
        if self.species == "Unknown":
            return "No species-specific markers detected"
        markers = [m["gene"] for m in self.detected_markers]
        return f"Identified as {self.species} based on: {', '.join(markers)}"


def _db_version() -> str:
    try:
        return hashlib.sha256(_MARKERS_FASTA.read_bytes()).hexdigest()[:8]
    except OSError:
        return "unknown"


def _load_rules() -> list[dict[str, Any]]:
    if not _MARKER_RULES.is_file():
        return []
    data = yaml.safe_load(_MARKER_RULES.read_text())
    return list(data.get("rules", []))


def _blast_contigs(contigs_fasta: str, threads: int = 4) -> list[dict[str, Any]]:
    from gside.config import which
    from gside.engine._env import require_bin
    from gside.engine.backends.blast import BlastBackend

    blastn = require_bin(
        "blastn",
        hint="blastn not found — install via pixi (conda) or point GSIDE_PIXI_BIN at a bin dir",
        resolver=which,
    )

    db = str(SPECIES_DB_DIR / "D1_marker" / "markers_blastdb")
    backend = BlastBackend(binary=blastn, threads=threads)
    engine_hits = backend.find(
        query=Path(contigs_fasta),
        db_path=db,
        evalue=1e-10,
        word_size=11,
    )
    hits = []
    for hit in engine_hits:
        seqid = hit.subject_id
        parts = seqid.split("~~~")
        gene = parts[1].lower() if len(parts) >= 2 else seqid.lower()
        hits.append(
            {
                "gene": gene,
                "identity": hit.identity,
                "coverage": round(hit.subject_coverage, 1),
                "ref_len": hit.subject_length,
                "seqid": seqid,
            }
        )
    return hits


def _best_per_gene(
    hits: list[dict[str, Any]],
    min_identity: float = 0.0,
    min_coverage: float = 0.0,
) -> dict[str, dict[str, Any]]:
    """Highest-identity hit per gene, optionally filtered by identity/coverage."""
    best: dict[str, dict[str, Any]] = {}
    for hit in hits:
        if hit["identity"] < min_identity or hit["coverage"] < min_coverage:
            continue
        gene = hit["gene"]
        if gene not in best or hit["identity"] > best[gene]["identity"]:
            best[gene] = hit
    return best


def _priority_species(rules: list[dict[str, Any]]) -> set[str]:
    """Species declared as priority targets by any rule."""
    return {rule["species"] for rule in rules if rule.get("priority_over")}


def _rule_match(
    rule: dict[str, Any], significant: dict[str, dict[str, Any]]
) -> tuple[list[str], list[str]]:
    """(matched genes passing the rule's min_identity, excluded genes present)."""
    genes = [g.lower() for g in rule.get("genes", [])]
    exclude = [g.lower() for g in rule.get("exclude_genes", [])]
    min_id = rule.get("min_identity", _HIGH_CONF)
    matched = [
        gene
        for gene in genes
        if gene in significant and significant[gene]["identity"] >= min_id and gene not in exclude
    ]
    excluded_present = [g for g in exclude if g in significant]
    return matched, excluded_present


def _pick_rule(
    rules: list[dict[str, Any]],
    significant: dict[str, dict[str, Any]],
    priority_species: set[str],
) -> tuple[str, str, str]:
    """(species, confidence, matched_rule) for the best satisfied rule."""
    best_species, best_conf, best_score, best_rule = "Unknown", "low", 0, ""
    for rule in rules:
        min_hits = rule.get("min_hits", 1)
        matched, excluded_present = _rule_match(rule, significant)
        if len(matched) < min_hits or excluded_present:
            continue
        score = len(matched)
        avg_id = sum(significant[g]["identity"] for g in matched) / len(matched)
        conf = "high" if avg_id >= _HIGH_CONF else "medium"
        rule_str = f"{'+'.join(matched)} ({len(matched)}/{min_hits} required)"

        is_priority = rule["species"] in priority_species
        if is_priority:
            # priority 物种：一旦命中即压制非 priority 物种（即使后者分数更高）
            if best_species not in priority_species or score > best_score:
                best_species, best_conf, best_score, best_rule = (
                    rule["species"],
                    conf,
                    score,
                    rule_str,
                )
        elif best_species not in priority_species and (
            score > best_score or (score == best_score and conf == "high")
        ):
            best_species, best_conf, best_score, best_rule = (
                rule["species"],
                conf,
                score,
                rule_str,
            )
    return best_species, best_conf, best_rule


def identify_multigene(contigs_fasta: str, threads: int = 4) -> MultiGeneResult:
    raw_hits = _blast_contigs(contigs_fasta, threads)
    rules = _load_rules()
    significant = _best_per_gene(raw_hits, _MIN_IDENTITY, _MIN_COVERAGE)

    result = MultiGeneResult(database_version=_db_version())
    result.all_hits = list(_best_per_gene(raw_hits).values())
    result.species, result.confidence, result.matched_rule = _pick_rule(
        rules, significant, _priority_species(rules)
    )

    if significant:
        result.detected_markers = [
            {"gene": g, "identity": h["identity"], "coverage": h["coverage"]}
            for g, h in sorted(significant.items(), key=lambda x: -x[1]["identity"])
        ]

    if not significant and raw_hits:
        near = [h for h in raw_hits if h["identity"] >= 75 and h["coverage"] >= 20]
        if near:
            result.notes.append(
                f"Near-threshold hits detected (identity 75-85%): "
                f"{', '.join(h['gene'] for h in near[:5])}. "
                "Possible divergence or contamination; ANI recheck advised."
            )

    return result


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Multi-gene species identifier")
    parser.add_argument("contigs")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    res = identify_multigene(args.contigs)
    if args.json:
        print(json.dumps(res.to_dict(), ensure_ascii=False, indent=2))
    else:
        d = res.to_dict()
        print(f"Species: {d['species']} ({d['confidence']})")
        print(f"Method: {d['method']}, Rule: {d['matched_rule']}")
        for m in d.get("detected_markers", []):
            print(f"  {m['gene']}: {m['identity']}% identity, {m['coverage']}% coverage")
        for n in d.get("notes", []):
            print(f"  NOTE: {n}")
        print(d["interpretation"])


if __name__ == "__main__":
    main()
