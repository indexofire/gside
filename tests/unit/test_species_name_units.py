"""Species label canonicalization (bacmap `Genus_species` contract)."""

from __future__ import annotations

import pytest

from gside.analysis.species_name import (
    SpeciesName,
    apply_species_name,
    canonical_species,
    display_name,
    parse_species_name,
)


class TestPanelAndMashNames:
    def test_space_binomial_with_strain(self):
        parsed = parse_species_name("Streptococcus agalactiae 2603V/R")
        assert parsed == SpeciesName(
            label="Streptococcus_agalactiae",
            display="Streptococcus agalactiae",
            strain="2603V/R",
        )

    def test_underscore_name_with_accession(self):
        parsed = parse_species_name("Streptococcus_agalactiae_GCF_949788565.1")
        assert parsed.label == "Streptococcus_agalactiae"
        assert parsed.display == "Streptococcus agalactiae"
        assert parsed.strain == "GCF_949788565.1"

    def test_space_binomial_without_strain(self):
        parsed = parse_species_name("Escherichia coli")
        assert parsed == SpeciesName(
            label="Escherichia_coli", display="Escherichia coli", strain=""
        )

    def test_ani_test_style_name(self):
        assert canonical_species("Escherichia coli K12") == "Escherichia_coli"

    def test_serogroup_suffix_is_carried_as_strain(self):
        parsed = parse_species_name("Vibrio cholerae O1")
        assert parsed.label == "Vibrio_cholerae"
        assert parsed.strain == "O1"


class TestMarkerRegistryLabels:
    @pytest.mark.parametrize(
        "label",
        [
            "Salmonella_enterica",
            "Vibrio_parahaemolyticus",
            "Escherichia_coli",
            "Shigella",
            "Shigella_flexneri",
            "Streptococcus_agalactiae",
        ],
    )
    def test_binomial_labels_are_idempotent(self, label):
        assert canonical_species(label) == label

    def test_registry_label_display_name(self):
        assert display_name("Vibrio_parahaemolyticus") == "Vibrio parahaemolyticus"
        assert display_name("") == "Unknown"

    def test_serotype_variant_label_is_not_reduced(self):
        # marker_rules.yaml emits `Vibrio_cholerae_O1`; bacmap accepts
        # multi-segment labels, so display expansion must not drop the suffix.
        assert display_name("Vibrio_cholerae_O1") == "Vibrio cholerae O1"
        assert display_name("Shigella_flexneri") == "Shigella flexneri"


class TestSourmashNames:
    def test_accession_prefix_is_stripped(self):
        parsed = parse_species_name("GCF_001457695.1 Campylobacter jejuni")
        assert parsed.label == "Campylobacter_jejuni"
        assert parsed.display == "Campylobacter jejuni"
        assert parsed.strain == ""

    def test_gca_accession_prefix_is_stripped(self):
        assert canonical_species("GCA_014305285.1 Campylobacter jejuni") == ("Campylobacter_jejuni")

    def test_bare_accession_is_unknown(self):
        assert canonical_species("GCF_001457695.1") == "Unknown"


class TestGtdbNames:
    def test_species_rank_prefix_stripped(self):
        assert canonical_species("s__Vibrio parahaemolyticus") == "Vibrio_parahaemolyticus"

    def test_full_lineage_uses_last_rank(self):
        lineage = (
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;"
            "f__Streptococcaceae;g__Streptococcus;s__Streptococcus agalactiae"
        )
        assert canonical_species(lineage) == "Streptococcus_agalactiae"

    def test_empty_species_rank_is_unknown(self):
        assert canonical_species("s__") == "Unknown"

    def test_lineage_without_species_rank_is_unknown(self):
        assert canonical_species("d__Bacteria;p__Firmicutes;g__Clostridium") == "Unknown"


class TestUnknownHandling:
    @pytest.mark.parametrize(
        "raw",
        ["", "  ", "Unknown", "unknown", "Mixed/Unknown", "n/a", "unidentified"],
    )
    def test_placeholders_become_unknown(self, raw):
        assert canonical_species(raw) == "Unknown"

    @pytest.mark.parametrize(
        "raw",
        ["Campylobacter sp. CCS1377", "Vibrio sp", "Listeria spp.", "Uncultured bacterium"],
    )
    def test_species_group_placeholders_become_unknown(self, raw):
        assert canonical_species(raw) == "Unknown"

    def test_unknown_has_no_strain(self):
        assert parse_species_name("Unknown").strain == ""

    @pytest.mark.parametrize("raw", ["_", ".", "/", " - "])
    def test_separator_only_is_unknown(self, raw):
        assert canonical_species(raw) == "Unknown"


class TestIdempotence:
    @pytest.mark.parametrize(
        "raw",
        [
            "Streptococcus agalactiae 2603V/R",
            "s__Vibrio parahaemolyticus",
            "Vibrio_parahaemolyticus",
            "Salmonella",
        ],
    )
    def test_canonical_is_idempotent(self, raw):
        once = canonical_species(raw)
        assert canonical_species(once) == once


class TestApplySpeciesName:
    def test_writes_label_raw_and_strain(self):
        result = apply_species_name({}, "Streptococcus agalactiae 2603V/R")
        assert result["species"] == "Streptococcus_agalactiae"
        assert result["species_raw"] == "Streptococcus agalactiae 2603V/R"
        assert result["strain"] == "2603V/R"

    def test_omits_raw_when_already_canonical(self):
        result = apply_species_name({}, "Streptococcus_agalactiae")
        assert result["species"] == "Streptococcus_agalactiae"
        assert "species_raw" not in result
        assert "strain" not in result

    def test_unknown_keeps_single_field(self):
        result = apply_species_name({"confidence": "low"}, "Unknown")
        assert result == {"confidence": "low", "species": "Unknown"}


class TestGtdbDisambiguationSuffix:
    def test_genus_suffix_is_dropped(self):
        parsed = parse_species_name("s__Campylobacter_D jejuni")
        assert parsed.label == "Campylobacter_jejuni"
        assert parsed.display == "Campylobacter jejuni"

    def test_species_suffix_is_dropped(self):
        assert canonical_species("Campylobacter jejuni_D") == "Campylobacter_jejuni"

    def test_abbreviated_genus_suffix(self):
        assert canonical_species("s__Vibrio_A cholerae") == "Vibrio_cholerae"

    def test_gtdb_placeholder_species_is_unknown(self):
        assert canonical_species("s__Dongia_A sp036004205") == "Unknown"

    def test_strain_tokens_are_preserved(self):
        assert parse_species_name("Vibrio cholerae O1").strain == "O1"
        assert parse_species_name("Escherichia coli K12").strain == "K12"
