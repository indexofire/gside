"""Shadow/parity checks for the optional gseek engine.

Run locally with the sibling checkout on PYTHONPATH, e.g.::

    PYTHONPATH=../gseek/src pixi run -e dev pytest tests/integration/test_gseek_parity.py -q

Without gseek installed the module is skipped, so the normal gside suite does
not gain an optional-install dependency.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("gseek")

from gseek.io import read_fasta  # noqa: E402

from gside.analysis.multigene_identifier import identify_multigene  # noqa: E402
from gside.analysis.read_screening import screen_reads  # noqa: E402
from gside.config import which  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
_MARKERS = _REPO / "data" / "db" / "D1_marker" / "markers.fasta"
_GENOME = _REPO / "data" / "db" / "D2_ani" / "genomes" / "GCF_000006745.1_ASM674v1_genomic.fna"

_READY = (
    _MARKERS.is_file()
    and _GENOME.is_file()
    and which("blastn") is not None
    and which("minimap2") is not None
)


@pytest.mark.skipif(not _READY, reason="marker DB or aligners missing")
class TestGseekParity:
    def test_marker_blast_matches_legacy_species_call(self, tmp_path, monkeypatch):
        monkeypatch.setenv("GSEEK_CACHE_DIR", str(tmp_path / "gseek-cache"))

        monkeypatch.setenv("GSIDE_ENGINE", "legacy")
        legacy = identify_multigene(str(_GENOME))
        monkeypatch.setenv("GSIDE_ENGINE", "gseek")
        new = identify_multigene(str(_GENOME))

        assert new.species == legacy.species
        assert new.confidence == legacy.confidence
        assert {marker["gene"] for marker in new.detected_markers} == {
            marker["gene"] for marker in legacy.detected_markers
        }

    def test_read_screening_matches_legacy_aggregate(self, tmp_path, monkeypatch):
        marker = read_fasta(_MARKERS)[0]
        sequence = marker.sequence[100:250]
        reads = tmp_path / "reads.fq"
        reads.write_text(
            f"@r1\n{sequence}\n+\n{'I' * len(sequence)}\n@r2\n{sequence[20:150]}\n+\n{'I' * 130}\n",
            encoding="utf-8",
        )
        monkeypatch.setenv("GSEEK_CACHE_DIR", str(tmp_path / "gseek-cache"))

        monkeypatch.setenv("GSIDE_ENGINE", "legacy")
        legacy = screen_reads(reads, _MARKERS, threads=2)
        monkeypatch.setenv("GSIDE_ENGINE", "gseek")
        new = screen_reads(reads, _MARKERS, threads=2)

        old_by_gene = {hit["gene"]: hit for hit in legacy}
        new_by_gene = {hit["gene"]: hit for hit in new}
        assert old_by_gene
        assert set(new_by_gene) == set(old_by_gene)
        for gene, old in old_by_gene.items():
            current = new_by_gene[gene]
            assert current["identity"] == pytest.approx(old["identity"], abs=1.0)
            assert current["coverage"] == pytest.approx(old["coverage"], abs=2.0)
