"""Custom marker database discovery/validation (no external tools)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from gside import marker_db
from gside.marker_db import MarkerDatabase


def _db(
    root: Path,
    *,
    fasta: str | None = None,
    rules: list[dict] | None = None,
    metadata: dict | None = None,
    with_rules: bool = True,
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "markers.fasta").write_text(
        fasta
        or (
            ">markers~~~geneA~~~refA marker A\nACGTACGT\n"
            ">markers~~~geneB~~~refB marker B\nTTTTCCCC\n"
        )
    )
    if with_rules:
        (root / "marker_rules.yaml").write_text(
            yaml.safe_dump(
                {
                    "rules": rules
                    if rules is not None
                    else [{"species": "Test_species", "genes": ["geneA"], "min_hits": 1}]
                }
            )
        )
    if metadata is not None:
        (root / "metadata.yaml").write_text(yaml.safe_dump(metadata))
    return root


def test_load_custom_database_metadata(tmp_path: Path) -> None:
    root = _db(
        tmp_path / "my_markers",
        metadata={"name": "my_markers", "description": "test DB", "version": "7"},
    )
    db = MarkerDatabase.from_path(root)
    assert db.custom is True
    assert db.name == "my_markers"
    assert db.description == "test DB"
    assert db.version == "7"
    assert db.blast_db == root / "markers_blastdb"
    assert db.to_dict()["schema"] == marker_db.MARKER_DB_SCHEMA


def test_load_parent_db_dir(tmp_path: Path) -> None:
    root = _db(tmp_path / "custom-root" / "D1_marker")
    db = MarkerDatabase.from_path(root.parent)
    assert db.root == root


def test_content_version_without_metadata(tmp_path: Path) -> None:
    db = MarkerDatabase.from_path(_db(tmp_path / "db"))
    assert len(db.version) == 8
    assert db.version != "unknown"


def test_validate_ok(tmp_path: Path) -> None:
    assert MarkerDatabase.from_path(_db(tmp_path / "db")).validate() == []


def test_validate_missing_files(tmp_path: Path) -> None:
    root = tmp_path / "db"
    root.mkdir()
    with pytest.raises(ValueError, match="markers.fasta"):
        MarkerDatabase.from_path(root)


def test_validate_bad_fasta_header_and_unknown_rule_gene(tmp_path: Path) -> None:
    root = _db(
        tmp_path / "db",
        fasta=">bad_header\nACGT\n>markers~~~geneA~~~refA\nACGT\n",
        rules=[{"species": "X", "genes": ["missing"], "min_hits": 1}],
    )
    errors = MarkerDatabase.from_path(root).validate()
    assert any("markers~~~<gene>" in err for err in errors)
    assert any("unknown gene" in err for err in errors)


def test_validate_duplicate_fasta_id(tmp_path: Path) -> None:
    root = _db(
        tmp_path / "db",
        fasta=">markers~~~geneA~~~refA\nACGT\n>markers~~~geneA~~~refA\nTTTT\n",
    )
    assert any("duplicate" in err for err in MarkerDatabase.from_path(root).validate())


def test_resolve_env_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root = _db(tmp_path / "envdb", metadata={"name": "envdb"})
    monkeypatch.setenv("GSIDE_MARKER_DB", str(root))
    db = marker_db.resolve_marker_db()
    assert db is not None and db.name == "envdb"


def test_default_marker_db_is_bundled() -> None:
    db = marker_db.default_marker_db()
    assert db.custom is False
    assert db.name == marker_db.DEFAULT_MARKER_NAME
    assert db.fasta.is_file()


def test_identify_multigene_uses_custom_db(monkeypatch, tmp_path: Path) -> None:
    from gside.analysis import multigene_identifier as mg

    db = MarkerDatabase.from_path(
        _db(
            tmp_path / "db",
            metadata={"name": "custom_markers", "version": "9"},
            rules=[{"species": "Test_species", "genes": ["geneA"], "min_hits": 1}],
        )
    )
    seen = {}

    def fake_blast(contigs, threads=4, marker_db=None):
        seen["db"] = marker_db
        return [{"gene": "genea", "identity": 99.0, "coverage": 100.0}]

    monkeypatch.setattr(mg, "_blast_contigs", fake_blast)
    result = mg.identify_multigene("sample.fna", marker_db=db)
    assert seen["db"] is db
    assert result.species == "Test_species"
    assert result.database_name == "custom_markers"
    assert result.to_dict()["database"] == {"name": "custom_markers", "version": "9"}
