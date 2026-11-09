"""CPIN extension end-to-end: every reference genome must call its species.

The CPIN reference genomes live in the gitignored ``D2_ani/genomes/``
directory (copied in by the CPIN build); each test skips when its genome is
absent, so CI stays green without the large fixtures.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gside.config import which as _which

_REPO = Path(__file__).resolve().parents[2]
_GENOMES = _REPO / "data" / "db" / "D2_ani" / "genomes"
_CPIN_DB = _REPO / "data" / "db" / "D1_marker" / "cpin"
_HAS_BLASTN = _which("blastn") is not None

_EXPECTED: dict[str, str] = {
    # CPIN bacteria (RefSeq reference strains used to extract the markers)
    "GCF_000742195.1": "Brucella",
    "GCF_000195955.2": "Mycobacterium_tuberculosis",
    "GCF_000008685.2": "Borrelia_burgdorferi",
    "GCF_000834295.1": "Yersinia_pseudotuberculosis",
    "GCF_011064845.1": "Citrobacter",
    "GCF_000013125.1": "Anaplasmataceae",
    "GCF_988247085.1": "Escherichia_albertii",
    "GCF_057076405.1": "Plesiomonas_shigelloides",
    "GCF_002796765.2": "Vibrio_fluvialis",
    "GCF_000833335.1": "Francisella_tularensis",
    "GCF_000011705.1": "Burkholderia_mallei",
    "GCF_000194745.1": "Burkholderia_gladioli",
    "GCF_000009805.1": "Nocardia",
    "GCF_058804955.1": "Bartonella",
    "GCF_044059125.1": "Rickettsiaceae",
    "GCF_000191925.1": "Chlamydia_psittaci",
    "GCF_000008745.1": "Chlamydia_pneumoniae",
    "GCF_000009065.1": "Yersinia_pestis",
    # CPIN fungi (ITS barcodes)
    "GCF_000149245.1": "Cryptococcus",
    "GCF_000002655.1": "Aspergillus",
    "GCF_001477535.1": "Pneumocystis",
    "GCF_000149335.2": "Coccidioides",
    "GCF_000150115.1": "Histoplasma",
    "GCF_000003525.1": "Blastomyces",
    "GCF_000150705.2": "Paracoccidioides",
    "GCF_046127475.1": "Talaromyces_marneffei",
    "GCF_000835555.1": "Rhinocladiella",
    "GCF_000835475.1": "Cladophialophora",
    "GCF_000182965.3": "Candida",
    "GCF_050370095.1": "Fusarium",
    "GCF_000961545.1": "Sporothrix",
    # pre-existing coverage regression (bundled panel genomes)
    "GCF_000005845.2": "Escherichia_coli",
    "GCF_000006745.1": "Vibrio_cholerae",
    "GCF_000006945.2": "Salmonella",
    "GCF_000009005.1": "Staphylococcus_aureus",
    "GCF_000009085.1": "Campylobacter_jejuni",
    "GCF_000465235.1": "Campylobacter_coli",
}


def _genome(accession: str) -> Path | None:
    matches = sorted(_GENOMES.glob(f"{accession}*.fna"))
    return matches[0] if matches else None


# These reference assemblies omit the rDNA tandem repeat, so the ITS barcode
# cannot self-validate against them (blastn finds zero ITS hits). The markers
# come from type/clinical strain records and are validated by the inter-genus
# cross matrix instead.
_RDNA_ABSENT_ASSEMBLIES = {
    "GCF_046127475.1",  # Talaromyces marneffei ATCC 18224 (2024 assembly)
    "GCF_000835555.1",  # Rhinocladiella mackenziei CBS 650.93
}


@pytest.mark.skipif(not _HAS_BLASTN, reason="blastn not available")
@pytest.mark.parametrize("accession,expected", sorted(_EXPECTED.items()))
def test_reference_genome_calls_expected_species(accession: str, expected: str):
    genome = _genome(accession)
    if genome is None:
        pytest.skip(f"{accession} not provisioned locally")
    if accession in _RDNA_ABSENT_ASSEMBLIES:
        pytest.skip(f"{accession} assembly lacks the rDNA cluster (ITS not testable)")
    from gside.analysis.multigene_identifier import identify_multigene

    cpin_species = expected not in {
        "Escherichia_coli",
        "Vibrio_cholerae",
        "Salmonella",
        "Staphylococcus_aureus",
        "Campylobacter_jejuni",
        "Campylobacter_coli",
    }
    if accession == "GCF_000009065.1":  # pestis needs default markers to exclude ptb
        marker_db = [str(_REPO / "data/db/D1_marker/default"), str(_CPIN_DB)]
    elif cpin_species:
        marker_db = str(_CPIN_DB)
    else:
        marker_db = None
    result = (
        identify_multigene(str(genome), marker_db=marker_db).to_dict()
        if marker_db
        else identify_multigene(str(genome)).to_dict()
    )
    assert result["species"] == expected, (
        f"{accession} called {result['species']} "
        f"(confidence={result['confidence']}, rule={result['matched_rule']!r})"
    )
