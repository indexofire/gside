"""Multi-gene species identifier using the expanded markers BLAST database.

Replaces the single-gene species_identifier for routine use. Reads
marker_rules.yaml for gene combination logic and scoring.
"""

from __future__ import annotations

import hashlib
import json
import statistics
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from gside.analysis.read_screening import is_fastq, screen_reads
from gside.config import SPECIES_DB_DIR
from gside.marker_db import MarkerDatabase, resolve_marker_db

_DEFAULT_MARKER_DIR = SPECIES_DB_DIR / "D1_marker" / "default"
_MARKERS_FASTA = _DEFAULT_MARKER_DIR / "markers.fasta"
_MARKER_RULES = _DEFAULT_MARKER_DIR / "marker_rules.yaml"
_MARKER_METADATA = _DEFAULT_MARKER_DIR / "metadata.yaml"

_MIN_IDENTITY = 85.0
_MIN_COVERAGE = 60.0
_HIGH_CONF = 90.0
_MIN_DEPTH = 5.0
_MIN_RELATIVE_DEPTH = 0.10
_MIN_MAPQ = 20
_CONTROL_MIN_IDENTITY = 80.0
_CONTROL_MIN_COVERAGE = 50.0


@dataclass
class MultiGeneResult:
    species: str = "Unknown"
    confidence: str = "low"
    detected_markers: list[dict[str, Any]] = field(default_factory=list)
    all_hits: list[dict[str, Any]] = field(default_factory=list)
    matched_rule: str = ""
    method: str = "multigene"
    database_version: str = "unknown"
    database_name: str = "species_markers_v2"
    notes: list[str] = field(default_factory=list)
    control_gene: str = ""
    control_depth: float | None = None
    control_source: str = ""
    depth_gate: str = ""
    databases: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "species": self.species,
            "confidence": self.confidence,
            "method": self.method,
            "database": {"name": self.database_name, "version": self.database_version},
            "detected_markers": self.detected_markers,
            "matched_rule": self.matched_rule,
            "control_gene": self.control_gene,
            "control_depth": self.control_depth,
            "control_source": self.control_source,
            "depth_gate": self.depth_gate,
            "notes": self.notes,
            "interpretation": self._interpret(),
        }
        if self.databases:
            payload["databases"] = self.databases
        return payload

    def _interpret(self) -> str:
        if self.species == "Unknown":
            return "No species-specific markers detected"
        markers = [m["gene"] for m in self.detected_markers]
        return f"Identified as {self.species} based on: {', '.join(markers)}"


def _coerce_marker_dbs(
    marker_db: str | Path | MarkerDatabase | Sequence[str | Path | MarkerDatabase] | None,
) -> list[MarkerDatabase]:
    """Normalize the ``marker_db`` argument to a list of databases.

    ``None`` (bundled default), a single path/database, and a sequence of
    them are all accepted; a one-element sequence collapses to a single run
    so its output stays identical to the single-database path.
    """
    if isinstance(marker_db, MarkerDatabase):
        return [marker_db]
    if marker_db is None or isinstance(marker_db, (str, Path)):
        if marker_db is None:
            resolved = resolve_marker_db(None)
            return [resolved] if resolved is not None else []
        return [MarkerDatabase.from_path(marker_db)]
    dbs: list[MarkerDatabase] = []
    for item in marker_db:
        if item is None:
            continue
        if isinstance(item, MarkerDatabase):
            dbs.append(item)
        else:
            dbs.append(MarkerDatabase.from_path(item))
    return dbs


def _db_version(marker_db: MarkerDatabase | None = None) -> str:
    if marker_db is not None:
        return marker_db.version
    try:
        return hashlib.sha256(_MARKERS_FASTA.read_bytes()).hexdigest()[:8]
    except OSError:
        return "unknown"


def _load_rules(marker_db: MarkerDatabase | None = None) -> list[dict[str, Any]]:
    path = marker_db.rules if marker_db is not None else _MARKER_RULES
    if not path.is_file():
        return []
    data = yaml.safe_load(path.read_text())
    return list(data.get("rules", []))


def _load_metadata(marker_db: MarkerDatabase | None = None) -> dict[str, Any]:
    path = marker_db.root / "metadata.yaml" if marker_db is not None else _MARKER_METADATA
    if not path.is_file():
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError:
        return {}
    return data if isinstance(data, dict) else {}


def _control_config(
    marker_db: MarkerDatabase | None = None,
) -> tuple[dict[str, list[str]], dict[str, Any]]:
    """(species -> control genes, reads-mode defaults) from metadata.yaml."""
    metadata = _load_metadata(marker_db)
    raw = metadata.get("control_genes") or {}
    controls: dict[str, list[str]] = {}
    if isinstance(raw, dict):
        for species, genes in raw.items():
            if isinstance(genes, str):
                genes = [genes]
            if isinstance(genes, list):
                controls[str(species)] = [str(g).lower() for g in genes if str(g).strip()]
    defaults = metadata.get("reads") if isinstance(metadata.get("reads"), dict) else {}
    return controls, dict(defaults or {})


@dataclass
class _DepthGate:
    hits: list[dict[str, Any]]
    control_depth: float | None = None
    control_gene: str = ""
    control_source: str = ""
    filtered: list[dict[str, Any]] = field(default_factory=list)


def _control_depth(
    hits: list[dict[str, Any]],
    control_genes: list[str],
) -> tuple[float | None, str, str]:
    """Control mean depth, label, and source (dedicated/fallback).

    A dedicated control is the median mean depth of the configured
    single-copy control gene(s).  If none is detected, fall back to the
    sample median of marker genes that pass the global identity/coverage
    floor; if that is empty too, the relative gate is disabled.
    """
    # Control alleles can be divergent (e.g. C. jejuni lpxA vs C. coli
    # lpxA_2).  For depth, use the deepest allele per gene; identity ranking
    # would pick a few exact reads from a minor contaminant instead of the
    # sample's real single-copy control.
    best_depth: dict[str, dict[str, Any]] = {}
    for hit in hits:
        depth = float(hit.get("mean_depth") or 0.0)
        if depth <= 0:
            continue
        gene = str(hit["gene"])
        current = best_depth.get(gene)
        if current is None or depth > float(current.get("mean_depth") or 0.0):
            best_depth[gene] = hit

    dedicated: list[tuple[str, float]] = []
    for gene in control_genes:
        hit = best_depth.get(gene)
        if hit is None:
            continue
        if hit["identity"] < _CONTROL_MIN_IDENTITY or hit["coverage"] < _CONTROL_MIN_COVERAGE:
            continue
        depth = float(hit.get("mean_depth") or 0.0)
        if depth > 0:
            dedicated.append((gene, depth))
    if dedicated:
        depths = [depth for _, depth in dedicated]
        label = "+".join(gene for gene, _ in dedicated)
        return statistics.median(depths), label, "dedicated"

    fallback = [
        float(hit["mean_depth"])
        for hit in best_depth.values()
        if float(hit.get("mean_depth") or 0.0) > 0
        and hit["identity"] >= _MIN_IDENTITY
        and hit["coverage"] >= _MIN_COVERAGE
    ]
    if fallback:
        return statistics.median(fallback), "sample_median", "fallback"
    return None, "", ""


def _apply_depth_gate(
    hits: list[dict[str, Any]],
    control_genes: list[str],
    min_depth: float,
    min_relative_depth: float,
) -> _DepthGate:
    """Gate reads-mode hits on absolute and control-relative mean depth."""
    gate = _DepthGate(hits=[])
    if not any("mean_depth" in hit for hit in hits):
        # Assembly (BLAST) hits have no depth; preserve legacy behaviour.
        gate.hits = list(hits)
        return gate

    control_depth, label, source = _control_depth(hits, control_genes)
    gate.control_depth = control_depth
    gate.control_gene = label
    gate.control_source = source
    if control_depth is None or control_depth <= 0:
        # No relative control: still enforce the absolute depth floor.
        for raw_hit in hits:
            hit = dict(raw_hit)
            depth = hit.get("mean_depth")
            if depth is None or float(depth) >= min_depth:
                gate.hits.append(hit)
            else:
                gate.filtered.append(hit)
        return gate

    for raw_hit in hits:
        hit = dict(raw_hit)
        depth = hit.get("mean_depth")
        if depth is None:
            gate.hits.append(hit)
            continue
        relative = float(depth) / control_depth
        hit["relative_depth_pct"] = round(100.0 * relative, 1)
        if float(depth) >= min_depth and relative >= min_relative_depth:
            gate.hits.append(hit)
        else:
            gate.filtered.append(hit)
    return gate


def _legacy_blast_hits(
    contigs_fasta: str,
    threads: int,
    marker_db: MarkerDatabase | None = None,
) -> list[Any]:
    from gside.config import which
    from gside.engine._env import require_bin
    from gside.engine.backends.blast import BlastBackend

    blastn = require_bin(
        "blastn",
        hint="blastn not found — install via pixi (conda) or point GSIDE_PIXI_BIN at a bin dir",
        resolver=which,
    )

    if marker_db is not None:
        db = str(marker_db.ensure_blast_db())
    else:
        db = str(SPECIES_DB_DIR / "D1_marker" / "default" / "markers_blastdb")
    backend = BlastBackend(binary=blastn, threads=threads)
    return backend.find(
        query=Path(contigs_fasta),
        db_path=db,
        evalue=1e-10,
        word_size=11,
    )


def _blast_contigs(
    contigs_fasta: str,
    threads: int = 4,
    marker_db: MarkerDatabase | None = None,
) -> list[dict[str, Any]]:
    from gside.engine.gseek_bridge import blast_hits, engine, log_blast_shadow

    markers_fasta = marker_db.fasta if marker_db is not None else _MARKERS_FASTA

    def legacy() -> list[Any]:
        return _legacy_blast_hits(contigs_fasta, threads, marker_db)

    def new() -> list[Any]:
        return blast_hits(
            Path(contigs_fasta),
            markers_fasta,
            evalue=1e-10,
            word_size=11,
            threads=threads,
        )

    def compare(old: list[Any], new_hits: list[Any]) -> bool:
        log_blast_shadow(old, new_hits)
        return True

    engine_hits = engine().run(legacy, new, compare=compare)

    hits = []
    for hit in engine_hits:
        seqid = hit.subject_id
        from gside.marker_db import gene_from_seqid

        gene = gene_from_seqid(seqid).lower()
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
        if gene not in best or (hit["identity"], hit["coverage"]) > (
            best[gene]["identity"],
            best[gene]["coverage"],
        ):
            best[gene] = hit
    return best


def _priority_species(rules: list[dict[str, Any]]) -> set[str]:
    """Species declared as priority targets by any rule."""
    return {rule["species"] for rule in rules if rule.get("priority_over")}


def _hits_by_gene(hits: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """All hits grouped by gene, so each rule can apply its own thresholds."""
    grouped: dict[str, list[dict[str, Any]]] = {}
    for hit in hits:
        grouped.setdefault(hit["gene"], []).append(hit)
    return grouped


def _rule_match(
    rule: dict[str, Any],
    by_gene: dict[str, list[dict[str, Any]]],
    significant: dict[str, dict[str, Any]],
) -> tuple[list[str], list[str], list[float]]:
    """(matched genes, excluded genes present, identity per matched gene).

    ``min_identity`` / ``min_coverage`` are per-rule (defaults: 90 / 60) and
    are applied to every raw hit of the gene, so a lower-identity hit with
    better coverage can still satisfy the rule.
    """
    min_id = rule.get("min_identity", _HIGH_CONF)
    min_cov = rule.get("min_coverage", _MIN_COVERAGE)
    exclude = [g.lower() for g in rule.get("exclude_genes", [])]

    matched: list[str] = []
    identities: list[float] = []
    for gene in (g.lower() for g in rule.get("genes", [])):
        if gene in exclude:
            continue
        passing = [
            hit
            for hit in by_gene.get(gene, [])
            if hit["identity"] >= min_id and hit["coverage"] >= min_cov
        ]
        if passing:
            matched.append(gene)
            identities.append(max(hit["identity"] for hit in passing))

    # Exclusions keep the global detection threshold (identity + coverage).
    excluded_present = [g for g in exclude if g in significant]
    return matched, excluded_present, identities


def _pick_rule(
    rules: list[dict[str, Any]],
    by_gene: dict[str, list[dict[str, Any]]],
    significant: dict[str, dict[str, Any]],
    priority_species: set[str],
) -> tuple[str, str, str]:
    """(species, confidence, matched_rule) for the best satisfied rule."""
    best_species, best_conf, best_score, best_rule = "Unknown", "low", 0, ""
    for rule in rules:
        min_hits = rule.get("min_hits", 1)
        matched, excluded_present, identities = _rule_match(rule, by_gene, significant)
        if len(matched) < min_hits or excluded_present:
            continue
        score = len(matched)
        avg_id = sum(identities) / len(identities)
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


def identify_multigene(
    contigs_fasta: str | list[str],
    threads: int = 4,
    marker_db: str | Path | MarkerDatabase | Sequence[str | Path | MarkerDatabase] | None = None,
    min_depth: float | None = None,
    min_relative_depth: float | None = None,
    min_mapq: int | None = None,
) -> MultiGeneResult:
    """Identify from an assembly or one/more FASTQ files (mates combined).

    Pass ``marker_db`` to use a custom marker database directory, or set
    ``GSIDE_MARKER_DB``. Without either, the bundled ``D1_marker`` database is
    used.

    Reads mode applies a gspec-style abundance gate: a marker must have
    ``mean_depth >= min_depth`` and be at least ``min_relative_depth`` of a
    conserved single-copy control gene.  Assembly mode keeps the legacy
    BLAST identity/coverage rules (depth is not available from BLAST hits).
    """
    dbs = _coerce_marker_dbs(marker_db)
    multi = len(dbs) > 1
    db = dbs[0] if len(dbs) == 1 else None
    files = [contigs_fasta] if isinstance(contigs_fasta, str) else list(contigs_fasta)
    reads_mode = is_fastq(files[0])
    reads: str | list[str] = files[0] if len(files) == 1 else files

    if multi:
        # Union semantics: rules from every database compete over the pooled
        # hits; control genes and reads defaults merge first-wins so the
        # selection order stays deterministic.
        rules: list[dict[str, Any]] = []
        controls: dict[str, list[str]] = {}
        read_defaults: dict[str, Any] = {}
        for entry in dbs:
            rules.extend(_load_rules(entry))
            entry_controls, entry_defaults = _control_config(entry)
            for species_key, genes in entry_controls.items():
                controls.setdefault(species_key, genes)
            for key, value in entry_defaults.items():
                read_defaults.setdefault(key, value)
    else:
        controls, read_defaults = _control_config(db)
        rules = _load_rules(db) if db is not None else _load_rules()

    if min_depth is None:
        min_depth = float(read_defaults.get("min_depth", _MIN_DEPTH))
    if min_relative_depth is None:
        min_relative_depth = float(read_defaults.get("min_relative_depth", _MIN_RELATIVE_DEPTH))
    if min_mapq is None:
        min_mapq = int(read_defaults.get("min_mapq", _MIN_MAPQ))

    raw_hits: list[dict[str, Any]] = []
    if multi:
        for entry in dbs:
            if reads_mode:
                raw_hits.extend(screen_reads(reads, entry.fasta, threads, min_mapq=min_mapq))
            else:
                raw_hits.extend(_blast_contigs(files[0], threads, entry))
    elif reads_mode:
        markers_fasta = db.fasta if db is not None else _MARKERS_FASTA
        raw_hits = screen_reads(reads, markers_fasta, threads, min_mapq=min_mapq)
    elif db is not None:
        raw_hits = _blast_contigs(files[0], threads, db)
    else:
        raw_hits = _blast_contigs(files[0], threads)

    # Pick the control gene set from the preliminary (ungated) call, then gate.
    prelim_significant = _best_per_gene(raw_hits, _MIN_IDENTITY, _MIN_COVERAGE)
    prelim_species, _, _ = _pick_rule(
        rules,
        _hits_by_gene(raw_hits),
        prelim_significant,
        _priority_species(rules),
    )
    control_genes = controls.get(prelim_species, [])
    if reads_mode:
        gate = _apply_depth_gate(raw_hits, control_genes, min_depth, min_relative_depth)
    else:
        gate = _DepthGate(hits=list(raw_hits))
    hits = gate.hits
    significant = _best_per_gene(hits, _MIN_IDENTITY, _MIN_COVERAGE)

    if multi:
        databases_meta = [
            {"name": entry.name, "version": _db_version(entry), "path": str(entry.root)}
            for entry in dbs
        ]
        database_name = "+".join(d["name"] for d in databases_meta)
        database_version = "+".join(d["version"] for d in databases_meta)
    else:
        databases_meta = []
        database_name = db.name if db is not None else "species_markers_v2"
        database_version = _db_version(db) if db is not None else _db_version()

    result = MultiGeneResult(
        database_version=database_version,
        database_name=database_name,
        databases=databases_meta,
        control_gene=gate.control_gene,
        control_depth=round(gate.control_depth, 2) if gate.control_depth is not None else None,
        control_source=gate.control_source,
        depth_gate=(
            (
                f"min_depth={min_depth:g}x; "
                f"min_relative_depth={min_relative_depth:.0%}; "
                f"min_mapq={min_mapq}; control={gate.control_gene or 'none'}"
            )
            if reads_mode
            else ""
        ),
    )
    result.all_hits = list(_best_per_gene(hits).values())
    result.species, result.confidence, result.matched_rule = _pick_rule(
        rules, _hits_by_gene(hits), significant, _priority_species(rules)
    )
    if result.species == "Shigella":
        result.notes.append(
            "Shigella is reported at genus level only; species/serotype "
            "discrimination is not performed here. Use gspec for Shigella "
            "serotyping/typing."
        )

    all_control_genes = {gene for genes in controls.values() for gene in genes}
    rule_genes = {str(gene).lower() for rule in rules for gene in rule.get("genes", [])}
    # A control that is also a rule gene (e.g. gdh for S. suis) must stay
    # visible: it is the species marker itself, not a pure depth control.
    display_control_genes = all_control_genes - rule_genes
    if significant:
        for gene, hit in sorted(
            significant.items(), key=lambda item: (-item[1]["identity"], item[0])
        ):
            if gene in display_control_genes:
                continue
            marker = {
                "gene": gene,
                "identity": hit["identity"],
                "coverage": hit["coverage"],
            }
            for key in ("mean_depth", "reads", "relative_depth_pct"):
                if key in hit:
                    marker[key] = hit[key]
            result.detected_markers.append(marker)

    if reads_mode:
        if gate.control_depth is None:
            result.notes.append(
                "No control gene detected; relative-depth gate disabled "
                f"(absolute min_depth={min_depth:g}x still applied)."
            )
        elif gate.control_source == "fallback":
            result.notes.append(
                f"No dedicated control for {prelim_species}; used sample median "
                f"depth ({gate.control_depth:.1f}x) as the relative-depth control."
            )
        if gate.filtered:
            shown = sorted(
                (h for h in gate.filtered if h["gene"] not in all_control_genes),
                key=lambda h: -(h.get("mean_depth") or 0),
            )[:8]
            if shown:
                result.notes.append(
                    "Depth-gated out: "
                    + ", ".join(
                        f"{h['gene']}({h.get('mean_depth', 0):.1f}x/"
                        f"{h.get('relative_depth_pct', 0):.1f}%)"
                        for h in shown
                    )
                    + (" ..." if len(gate.filtered) > len(shown) else "")
                )

    if not significant and hits:
        near = sorted(
            (h for h in hits if h["identity"] >= 75 and h["coverage"] >= 20),
            key=lambda h: h["gene"],
        )
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
