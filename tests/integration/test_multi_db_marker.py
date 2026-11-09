"""Multi-database selection: real blastn across two custom marker databases.

Integration tests on purpose — they need the conda-provisioned ``blastn``
(resolved via ``GSIDE_PIXI_BIN`` or the pixi env). Resolution logic itself is
covered hermetically in ``tests/unit/test_db_registry_units.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gside.config import which as _which

_REPO = Path(__file__).resolve().parents[2]
CONTIGS = _REPO / "data" / "db" / "D1_marker" / "markers.fasta"
_HAS_BLASTN = _which("blastn") is not None


def _gene_seq(prefix: str) -> str:
    blocks = CONTIGS.read_text().split(">")
    block = next(b for b in blocks if b.startswith(prefix))
    return block.split("\n", 1)[1].replace("\n", "")


def _custom_db(root: Path, name: str, gene: str, species: str, seq: str) -> Path:
    path = root / "marker" / name
    path.mkdir(parents=True)
    (path / "markers.fasta").write_text(f">markers~~~{gene}~~~x marker {gene}\n{seq}\n")
    (path / "marker_rules.yaml").write_text(
        f"rules:\n- species: {species}\n  genes: [{gene}]\n  min_hits: 1\n"
    )
    return path


@pytest.mark.skipif(not _HAS_BLASTN, reason="blastn not available")
class TestMultiDatabaseMarker:
    def test_two_dbs_call_species_across_the_union(self, tmp_path):
        """Each database defines one species; the union run reports both genes
        and arbitrates to the higher-identity rule."""
        from gside.analysis.multigene_identifier import identify_multigene

        db_a = _custom_db(
            tmp_path, "dbA", "inva", "Salmonella_enterica", _gene_seq("markers~~~inva~~~")
        )
        db_b = _custom_db(
            tmp_path, "dbB", "tlh", "Vibrio_parahaemolyticus", _gene_seq("markers~~~tlh~~~")
        )
        q = tmp_path / "union.fna"
        q.write_text(
            f">inva_q\n{_gene_seq('markers~~~inva~~~')}\n>tlh_q\n{_gene_seq('markers~~~tlh~~~')}\n"
        )
        payload = identify_multigene(str(q), marker_db=[str(db_a), str(db_b)]).to_dict()
        assert payload["databases"] == [
            {"name": "dbA", "version": payload["databases"][0]["version"], "path": str(db_a)},
            {"name": "dbB", "version": payload["databases"][1]["version"], "path": str(db_b)},
        ]
        assert payload["database"]["name"] == "dbA+dbB"
        # Both genes are detected; the species is one of the two rules
        # (arbitration over identical identities keeps the first-seen rule).
        assert {m["gene"] for m in payload["detected_markers"]} == {"inva", "tlh"}
        assert payload["species"] in {"Salmonella_enterica", "Vibrio_parahaemolyticus"}

    def test_cli_d_flag_selects_named_databases(self, tmp_path):
        """`-d` names resolve under the family dir of the -D root."""
        from click.testing import CliRunner

        from gside.cli import main

        _custom_db(tmp_path, "food", "inva", "Salmonella_enterica", _gene_seq("markers~~~inva~~~"))
        q = tmp_path / "inva.fna"
        q.write_text(f">inva_q\n{_gene_seq('markers~~~inva~~~')}\n")
        result = CliRunner().invoke(
            main,
            ["scan", str(q), "-D", str(tmp_path), "-d", "food", "-f", "json"],
        )
        assert result.exit_code == 0, result.output
        assert '"species": "Salmonella_enterica"' in result.output
        assert '"name": "food"' in result.output
