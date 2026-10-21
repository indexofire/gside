"""Pure unit tests for CLI rendering, arbitration and the marker table shape."""

from __future__ import annotations

from gside.cli import TABLE_COLUMNS, _arbitrate, _render_table


def _payload() -> dict:
    return {
        "contigs": "a.fna",
        "methods": {
            "marker": {
                "species": "Salmonella",
                "confidence": "high",
                "matched_rule": "inva (1/1 required)",
            }
        },
        "verdict": {"species": "Salmonella", "confidence": "high", "basis": ["marker"]},
    }


class TestRenderTable:
    def test_tsv_single_row(self):
        out = _render_table([_payload()], markdown=False).splitlines()
        assert out[0].split("\t") == TABLE_COLUMNS
        row = out[1].split("\t")
        assert row[0] == "a.fna" and row[1] == "Salmonella" and row[6] == "inva (1/1 required)"

    def test_md_shape(self):
        out = _render_table([_payload()], markdown=True).splitlines()
        assert out[0].startswith("| contigs |") and out[1].count("---") == 14
        assert "Salmonella" in out[2]

    def test_error_cell(self):
        payload = {
            "contigs": "b.fna",
            "methods": {"panel": {"method": "panel", "error": "db missing"}},
            "verdict": {"species": "Unknown", "confidence": "low", "basis": ["panel"]},
        }
        row = _render_table([payload], markdown=False).splitlines()[1].split("\t")
        assert row[1] == "Unknown" and "db missing" in row[13]


class TestArbitrate:
    def test_ani_beats_marker(self):
        verdict = _arbitrate(
            {
                "marker": {"species": "X_marker", "confidence": "high"},
                "panel": {"result": {"species": "Y_ani", "confidence": "high"}},
            }
        )
        assert verdict["species"] == "Y_ani"
        assert verdict["basis"] == ["panel"]

    def test_errors_skipped(self):
        verdict = _arbitrate(
            {"panel": {"error": "db missing"}, "marker": {"species": "M", "confidence": "high"}}
        )
        assert verdict["species"] == "M"


class TestArbitrateRanking:
    def test_medium_beats_low_within_same_layer(self):
        verdict = _arbitrate(
            {
                "mash_refseq": {"result": {"species": "Mash_low", "confidence": "low"}},
                "panel": {"result": {"species": "Panel_medium", "confidence": "medium"}},
            }
        )
        assert verdict["species"] == "Panel_medium"

    def test_supporting_lists_agreeing_methods(self):
        verdict = _arbitrate(
            {
                "marker": {"species": "Streptococcus_agalactiae", "confidence": "high"},
                "panel": {
                    "result": {
                        "species": "Streptococcus_agalactiae",
                        "confidence": "high",
                        "strain": "2603V/R",
                    }
                },
                "mash_refseq": {
                    "result": {"species": "Staphylococcus_aureus", "confidence": "medium"}
                },
            }
        )
        assert verdict["species"] == "Streptococcus_agalactiae"
        assert verdict["basis"] == ["panel"]
        assert verdict["supporting"] == ["marker"]
        assert verdict["agreement"] == "2/3"
        assert verdict["strain"] == "2603V/R"

    def test_unknown_is_not_support(self):
        verdict = _arbitrate(
            {
                "marker": {"species": "Salmonella", "confidence": "high"},
                "panel": {"result": {"species": "Unknown", "confidence": "low"}},
            }
        )
        assert verdict["supporting"] == []
        assert verdict["agreement"] == "1/2"

    def test_display_name_expands_serotype_suffix(self):
        verdict = _arbitrate({"marker": {"species": "Vibrio_cholerae_O1", "confidence": "high"}})
        assert verdict["display_name"] == "Vibrio cholerae O1"

    def test_low_confidence_still_reported(self):
        verdict = _arbitrate({"marker": {"species": "Salmonella", "confidence": "low"}})
        assert (verdict["species"], verdict["confidence"]) == ("Salmonella", "low")
        assert verdict["agreement"] == "1/1"


class TestSingleModeVerdict:
    def test_carries_canonical_fields(self):
        from gside.cli import _verdict_for

        verdict = _verdict_for(
            {
                "panel": {
                    "result": {
                        "species": "Escherichia_coli",
                        "confidence": "high",
                        "strain": "K12",
                    }
                }
            }
        )
        assert verdict["species"] == "Escherichia_coli"
        assert verdict["display_name"] == "Escherichia coli"
        assert verdict["strain"] == "K12"
        assert verdict["supporting"] == []
        assert verdict["agreement"] == "1/1"


class TestCsvRendering:
    def test_species_csv_header_and_row(self, capsys):
        from gside.cli import _emit_species

        assert _emit_species(iter([_payload()]), 1, "csv", None) == 0
        lines = capsys.readouterr().out.splitlines()
        assert lines[0].startswith("contigs,verdict_species,")
        assert lines[1].startswith("a.fna,Salmonella,high,marker,")

    def test_species_csv_no_header(self, capsys):
        from gside.cli import _emit_species

        assert _emit_species(iter([_payload()]), 1, "csv", None, no_header=True) == 0
        lines = capsys.readouterr().out.splitlines()
        assert len(lines) == 1
        assert lines[0].startswith("a.fna,")
