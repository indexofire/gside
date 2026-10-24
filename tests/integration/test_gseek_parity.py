"""Shadow/parity checks for the optional gseek engine.

Run locally with the sibling checkout on PYTHONPATH, e.g.::

    PYTHONPATH=../gseek/src pixi run -e dev pytest tests/integration/test_gseek_parity.py -q

Without gseek installed the module is skipped, so the normal gside suite does
not gain an optional-install dependency.
"""

from __future__ import annotations

import random
import subprocess
from pathlib import Path

import pytest

pytest.importorskip("gseek")

from gseek.io import read_fasta  # noqa: E402

from gside.analysis.ani_identifier import identify_by_ani  # noqa: E402
from gside.analysis.multigene_identifier import identify_multigene  # noqa: E402
from gside.analysis.read_screening import screen_reads  # noqa: E402
from gside.config import which  # noqa: E402
from gside.engine import gseek_bridge  # noqa: E402

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


_MASH_DB = _REPO / "data" / "db" / "D3_mash" / "mash.msh"
_K12 = _REPO / "data" / "db" / "D2_ani" / "genomes" / "GCF_000005845.2_ASM584v2_genomic.fna"


@pytest.mark.skipif(
    not (_MASH_DB.is_file() and _K12.is_file()) or which("mash") is None,
    reason="mash DB or binary missing",
)
def test_mash_refseq_matches_legacy(tmp_path, monkeypatch):
    from gside.analysis.ani_identifier import identify_by_ani

    monkeypatch.setenv("GSEEK_CACHE_DIR", str(tmp_path / "gseek-cache"))
    monkeypatch.setenv("GSIDE_ENGINE", "legacy")
    legacy = identify_by_ani(str(_K12), mode="mash_refseq").result
    monkeypatch.setenv("GSIDE_ENGINE", "gseek")
    new = identify_by_ani(str(_K12), mode="mash_refseq").result

    assert new["species"] == legacy["species"]
    assert new["confidence"] == legacy["confidence"]
    assert new["top_hits"][0]["genome"] == legacy["top_hits"][0]["genome"]
    assert new["identity"] == pytest.approx(legacy["identity"], abs=1e-3)


@pytest.mark.skipif(which("sourmash") is None, reason="sourmash not available")
def test_sourmash_external_signature_gather(tmp_path, monkeypatch):
    random.seed(20261032)
    sequence = "".join(random.choice("ACGT") for _ in range(200_000))
    target = tmp_path / "target.fa"
    target.write_text(f">target\n{sequence}\n", encoding="utf-8")
    query = tmp_path / "query.fa"
    query.write_text(f">query\n{sequence[500:95_500]}\n", encoding="utf-8")
    query2 = tmp_path / "query2.fa"
    query2.write_text(f">query2\n{sequence[100_000:195_000]}\n", encoding="utf-8")
    reference = tmp_path / "refs.sig"
    subprocess.run(
        [
            "sourmash",
            "sketch",
            "dna",
            "-p",
            "k=31,scaled=1000",
            "-o",
            str(reference),
            str(target),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    monkeypatch.setenv("GSEEK_CACHE_DIR", str(tmp_path / "gseek-cache"))
    rows = gseek_bridge.sourmash_gather([query, query2], reference)
    assert rows
    assert 0.0 <= float(rows[0]["f_unique_weighted"]) <= 1.0


@pytest.mark.skipif(
    not (_MASH_DB.is_file() and _K12.is_file()) or which("mash") is None,
    reason="mash DB or binary missing",
)
def test_mash_reads_matches_legacy(tmp_path, monkeypatch):
    from gside.analysis.ani_identifier import identify_by_ani

    genome = "".join(
        line.strip() for line in _K12.read_text().splitlines() if not line.startswith(">")
    )
    rng = random.Random(20261033)
    reads = tmp_path / "reads.fq"
    with reads.open("w", encoding="utf-8") as handle:
        for index in range(5000):
            start = rng.randrange(0, len(genome) - 150)
            sequence = genome[start : start + 150]
            handle.write(f"@r{index}\n{sequence}\n+\n{'I' * 150}\n")

    monkeypatch.setenv("GSEEK_CACHE_DIR", str(tmp_path / "gseek-cache"))
    monkeypatch.setenv("GSIDE_ENGINE", "legacy")
    legacy = identify_by_ani(str(reads), mode="mash_refseq").result
    monkeypatch.setenv("GSIDE_ENGINE", "gseek")
    new = identify_by_ani(str(reads), mode="mash_refseq").result

    assert new["species"] == legacy["species"]
    assert new["confidence"] == legacy["confidence"]
    assert new["top_hits"][0]["genome"] == legacy["top_hits"][0]["genome"]
    assert new["identity"] == pytest.approx(legacy["identity"], abs=1e-3)


_PANEL_DB = _REPO / "data" / "db" / "D2_ani" / "panel.sketch"


@pytest.mark.skipif(
    not _PANEL_DB.is_dir() or which("skani") is None,
    reason="panel sketched DB or skani missing",
)
def test_panel_matches_legacy(tmp_path, monkeypatch):
    monkeypatch.setenv("GSEEK_CACHE_DIR", str(tmp_path / "gseek-cache"))
    monkeypatch.setenv("GSIDE_ENGINE", "legacy")
    legacy = identify_by_ani(str(_K12), mode="panel").result
    monkeypatch.setenv("GSIDE_ENGINE", "gseek")
    new = identify_by_ani(str(_K12), mode="panel").result
    assert new == legacy
