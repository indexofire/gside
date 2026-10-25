"""Marker database discovery, validation, and BLAST index management.

A marker database is a directory with:

``markers.fasta``
    Marker sequences. Headers must use ``markers~~~<gene>~~~<accession>``
    (the existing marker-header contract).
``marker_rules.yaml``
    Rules with the same schema as the bundled database.
``metadata.yaml`` (optional)
    ``name``, ``description``, and ``version`` for result provenance.

The bundled database lives at ``data/db/D1_marker``. Custom databases can be
any directory passed with ``--marker-db`` or ``$GSIDE_MARKER_DB``.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

MARKER_DB_SCHEMA = "gside.marker-db/1"
DEFAULT_MARKER_NAME = "species_markers_v2"
_REQUIRED_HEADER_PARTS = 3


@dataclass(frozen=True)
class MarkerDatabase:
    """A validated marker database on disk."""

    root: Path
    fasta: Path
    rules: Path
    name: str = DEFAULT_MARKER_NAME
    description: str = ""
    version: str = "unknown"
    metadata: Mapping[str, Any] = field(default_factory=dict)
    custom: bool = False

    @property
    def blast_db(self) -> Path:
        """BLAST database prefix (``makeblastdb -out`` value)."""
        configured = self.metadata.get("blast_db")
        return self.root / str(configured) if configured else self.root / "markers_blastdb"

    @classmethod
    def from_path(cls, path: str | Path) -> MarkerDatabase:
        """Load a marker DB from *path*.

        A parent directory containing exactly one ``D1_marker`` child is also
        accepted, so ``--db-dir`` roots work without an extra flag.
        """
        root = Path(path).expanduser().resolve()
        if not root.is_dir():
            raise ValueError(f"marker database path is not a directory: {root}")
        if not (root / "markers.fasta").is_file() and (root / "D1_marker").is_dir():
            root = root / "D1_marker"

        fasta = root / "markers.fasta"
        rules = root / "marker_rules.yaml"
        if not fasta.is_file() or not rules.is_file():
            raise ValueError(
                f"marker database {root} must contain markers.fasta and marker_rules.yaml"
            )

        metadata_path = root / "metadata.yaml"
        metadata: dict[str, Any] = {}
        if metadata_path.is_file():
            try:
                loaded = yaml.safe_load(metadata_path.read_text(encoding="utf-8")) or {}
            except yaml.YAMLError as exc:
                raise ValueError(f"invalid metadata.yaml: {exc}") from exc
            if not isinstance(loaded, dict):
                raise ValueError("metadata.yaml must contain a mapping")
            metadata = loaded

        schema = metadata.get("schema")
        if schema is not None and str(schema) != MARKER_DB_SCHEMA:
            raise ValueError(
                f"unsupported marker database schema {schema!r}; expected {MARKER_DB_SCHEMA}"
            )
        name = str(metadata.get("name") or root.name or DEFAULT_MARKER_NAME)
        description = str(metadata.get("description") or "")
        version = str(metadata.get("version") or _content_version(fasta, rules))
        return cls(
            root=root,
            fasta=fasta,
            rules=rules,
            name=name,
            description=description,
            version=version,
            metadata=metadata,
            custom=True,
        )

    def validate(self) -> list[str]:
        """Return human-readable validation errors; empty means valid."""
        errors: list[str] = []
        for label, path in (("markers.fasta", self.fasta), ("marker_rules.yaml", self.rules)):
            if not path.is_file():
                errors.append(f"missing {label}: {path}")

        genes: set[str] = set()
        if self.fasta.is_file():
            seen: set[str] = set()
            for line in self.fasta.read_text(encoding="utf-8", errors="replace").splitlines():
                if not line.startswith(">"):
                    continue
                seq_id = line[1:].split()[0] if line[1:].split() else ""
                if not seq_id:
                    errors.append("empty FASTA header")
                    continue
                if seq_id in seen:
                    errors.append(f"duplicate FASTA id: {seq_id}")
                seen.add(seq_id)
                parts = seq_id.split("~~~")
                if len(parts) < _REQUIRED_HEADER_PARTS or not parts[1]:
                    errors.append(f"FASTA id {seq_id!r} must be markers~~~<gene>~~~<accession>")
                else:
                    genes.add(parts[1].lower())

        rules: list[Any] = []
        if self.rules.is_file():
            try:
                data = yaml.safe_load(self.rules.read_text(encoding="utf-8")) or {}
            except yaml.YAMLError as exc:
                errors.append(f"invalid marker_rules.yaml: {exc}")
                data = {}
            if not isinstance(data, dict):
                errors.append("marker_rules.yaml must contain a mapping with a 'rules' list")
            else:
                rules = data.get("rules", [])
                if not isinstance(rules, list):
                    errors.append("marker_rules.yaml 'rules' must be a list")
                    rules = []
                for index, rule in enumerate(rules):
                    prefix = f"rules[{index}]"
                    if not isinstance(rule, dict):
                        errors.append(f"{prefix} must be a mapping")
                        continue
                    species = rule.get("species")
                    rule_genes = rule.get("genes")
                    if not isinstance(species, str) or not species:
                        errors.append(f"{prefix}.species must be a non-empty string")
                    if not isinstance(rule_genes, list) or not rule_genes:
                        errors.append(f"{prefix}.genes must be a non-empty list")
                        continue
                    for gene in rule_genes:
                        if not isinstance(gene, str) or not gene:
                            errors.append(f"{prefix}.genes contains an invalid name")
                        elif gene.lower() not in genes:
                            errors.append(f"{prefix} references unknown gene: {gene}")

        if not errors and not genes:
            errors.append("markers.fasta contains no marker sequences")
        if not errors and not rules:
            errors.append("marker_rules.yaml contains no rules")
        return errors

    def ensure_blast_db(self) -> Path:
        """Return the BLAST prefix, building it in-place on first use."""
        prefix = self.blast_db
        if Path(f"{prefix}.nhr").is_file() or Path(f"{prefix}.nin").is_file():
            return prefix
        if not self.root.is_dir() or not os.access(self.root, os.W_OK):
            raise RuntimeError(
                f"marker BLAST database is missing and {self.root} is not writable; "
                f"build it with: makeblastdb -in {self.fasta} -dbtype nucl -out {prefix}"
            )
        from gside.engine.backends.blast import BlastBackend

        BlastBackend(binary="blastn").make_db(self.fasta, prefix)
        return prefix

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": MARKER_DB_SCHEMA,
            "name": self.name,
            "description": self.description,
            "version": self.version,
            "path": str(self.root),
            "fasta": str(self.fasta),
            "rules": str(self.rules),
            "blast_db": str(self.blast_db),
            "custom": self.custom,
        }


def _content_version(fasta: Path, rules: Path) -> str:
    digest = hashlib.sha256()
    for path in (fasta, rules):
        try:
            digest.update(path.read_bytes())
        except OSError:
            return "unknown"
    return digest.hexdigest()[:8]


def resolve_marker_db(path: str | Path | None = None) -> MarkerDatabase | None:
    """Resolve an explicit path or ``$GSIDE_MARKER_DB``; None means bundled."""
    value = path or os.environ.get("GSIDE_MARKER_DB", "").strip() or None
    if value is None:
        return None
    return MarkerDatabase.from_path(value)


def default_marker_db(path: str | Path | None = None) -> MarkerDatabase:
    """Return an explicit/custom DB or the bundled ``D1_marker``."""
    resolved = resolve_marker_db(path)
    if resolved is not None:
        return resolved
    from gside.config import SPECIES_DB_DIR

    db = MarkerDatabase.from_path(SPECIES_DB_DIR / "D1_marker")
    return MarkerDatabase(
        root=db.root,
        fasta=db.fasta,
        rules=db.rules,
        name=DEFAULT_MARKER_NAME,
        description="Bundled gside marker database",
        version=db.version,
        metadata=db.metadata,
        custom=False,
    )


__all__ = [
    "MARKER_DB_SCHEMA",
    "MarkerDatabase",
    "default_marker_db",
    "resolve_marker_db",
]
