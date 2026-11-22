"""Salmonella species-level discrimination (S. enterica vs S. bongori).

Integration test over real RefSeq genomes placed in the (gitignored)
``D2_ani/genomes`` panel inputs: the two-species marker design —
``Salmonella_enterica`` (inva >= 93 identity band), ``Salmonella_bongori``
(inva + one bongori-specific marker, priority), and the bare ``Salmonella``
genus fallback (inva 90-93 without bongori markers). The empirical basis is
the 26-genome screen recorded in marker_rules.yaml notes (Fookes 2011 PLoS
Pathog e1002191 for the sbo/SPI-22 background). Skips when the verification
genomes are not provisioned (CI).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gside.analysis.multigene_identifier import identify_multigene

_GENOMES_DIR = Path(__file__).resolve().parents[2] / "data" / "db" / "D2_ani" / "genomes"

# 前缀 -> (期望物种, 期望 scope 说明)
_EXPECTED: dict[str, tuple[str, ...]] = {
    # S. enterica subsp. enterica serovar Typhimurium LT2 (inva ~99.8)
    "GCF_000006945.2": ("Salmonella_enterica",),
    # S. bongori NCTC 12419 reference (inva ~91.8 + asst/zadh ~100)
    "GCF_000252995.1": ("Salmonella_bongori",),
    # S. bongori field isolate (second lineage)
    "GCF_002035285.1": ("Salmonella_bongori",),
}


def _genome(accession: str):
    matches = sorted(_GENOMES_DIR.glob(f"{accession}*.fna"))
    return matches[0] if matches else None


def _markers(payload: dict) -> dict[str, float]:
    res = payload.get("result", payload)
    return {m["gene"]: m["identity"] for m in res.get("detected_markers", [])}


class TestSalmonellaSpeciesCalls:
    @pytest.mark.parametrize(("accession", "expected"), _EXPECTED.items())
    def test_species_call(self, accession: str, expected: tuple[str, ...]):
        genome = _genome(accession)
        if genome is None:
            pytest.skip(f"{accession} not provisioned under {_GENOMES_DIR}")
        payload = identify_multigene(str(genome), threads=4).to_dict()
        assert payload["species"] == expected[0], _markers(payload)

    def test_bongori_call_combines_genus_and_specific_markers(self):
        genome = _genome("GCF_000252995.1")
        if genome is None:
            pytest.skip("bongori reference genome not provisioned")
        payload = identify_multigene(str(genome), threads=4).to_dict()
        markers = _markers(payload)
        assert markers["inva"] >= 90.0, "genus signal must fire for the bongori rule"
        # enterica 规则被 exclude_genes 压制：即使 inva 过 90 也必须判 bongori
        assert markers.get("asst", 0) >= 90.0 or markers.get("zadh", 0) >= 90.0

    def test_enterica_rule_ships_with_bongori_rules(self):
        from gside.analysis.multigene_identifier import _load_rules

        rules = {rule["species"]: rule for rule in _load_rules()}
        assert set(rules) >= {"Salmonella", "Salmonella_enterica", "Salmonella_bongori"}
        enterica = rules["Salmonella_enterica"]
        assert enterica["min_identity"] == 90
        assert enterica["min_hits"] == 2
        assert enterica["require"] == ["inva"]
        assert sorted(enterica["exclude_genes"]) == ["asst", "zadh"]
        assert {"inva", "ssav", "ssec"} <= set(enterica["genes"])
        bongori = rules["Salmonella_bongori"]
        assert bongori["require"] == ["inva"]
        assert bongori["priority_over"] == [
            "Salmonella",
            "Salmonella_enterica",
        ]
        assert {"inva", "asst", "zadh"} <= set(bongori["genes"])
