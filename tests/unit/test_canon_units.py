"""Pure unit tests for species_canon (no importers in-repo, keep healthy)."""

from __future__ import annotations

from gside.analysis.species_canon import (
    ECOLI,
    SALMONELLA,
    SHIGELLA,
    UNKNOWN,
    VPARA,
    canonical_from_verdict,
    collapse_binomial,
)


class TestCanonicalFromVerdict:
    def test_empty(self):
        assert canonical_from_verdict("") == UNKNOWN
        assert canonical_from_verdict("   ") == UNKNOWN

    def test_not_prefix_splits_on_ipah(self):
        assert canonical_from_verdict("not_Shigella", "ipaH positive") == SHIGELLA
        assert canonical_from_verdict("not_Shigella", "negative") == ECOLI
        assert canonical_from_verdict("not_X") == ECOLI

    def test_needles(self):
        assert canonical_from_verdict("Salmonella enterica") == SALMONELLA
        assert canonical_from_verdict("Vibrio parahaemolyticus RIMD") == VPARA
        assert canonical_from_verdict("Shigella flexneri") == SHIGELLA
        assert canonical_from_verdict("E. coli K12") == ECOLI
        assert canonical_from_verdict("E.coli O157") == ECOLI
        assert canonical_from_verdict("DEC pathotype") == ECOLI

    def test_unknown_fallback(self):
        assert canonical_from_verdict("Listeria monocytogenes") == UNKNOWN

    def test_constants(self):
        assert (SALMONELLA, ECOLI, SHIGELLA, VPARA, UNKNOWN) == (
            "Salmonella",
            "E. coli",
            "Shigella",
            "V. parahaemolyticus",
            "unknown",
        )


class TestParseDbHeader:
    def test_all_arity_branches(self):
        from gside.utils import parse_db_header

        assert parse_db_header("db~~~g~~~acc~~~prod") == ("g", "acc", "prod", "")
        assert parse_db_header("db~~~g~~~acc") == ("g", "acc", "", "")
        assert parse_db_header("db~~~g") == ("g", "", "", "")
        assert parse_db_header("bare") == ("bare", "", "", "")


class TestCollapseBinomial:
    def test_not_prefix_passthrough(self):
        assert collapse_binomial("not_Shigella") == "not_Shigella"

    def test_mapping_and_fallback(self):
        assert collapse_binomial("Vibrio parahaemolyticus") == VPARA
        assert collapse_binomial("DEC") == ECOLI
        assert collapse_binomial("Unmapped thing") == "Unmapped thing"
