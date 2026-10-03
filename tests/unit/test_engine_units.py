"""Pure unit tests for engine/hits.py and engine/utils.py (no binaries)."""

from __future__ import annotations

import pytest

from gside.engine.hits import Hit
from gside.engine.utils import classify_allele, confidence_tier, merge_intervals


def _blast_line(**kw) -> str:
    cols = {
        "qseqid": "q1",
        "sseqid": "markers~~~inva~~~X",
        "pident": "99.5",
        "length": "400",
        "mismatch": "2",
        "gapopen": "0",
        "qstart": "10",
        "qend": "409",
        "sstart": "1",
        "send": "400",
        "evalue": "1e-100",
        "bitscore": "700",
        "qlen": "500",
        "slen": "450",
    }
    cols.update(kw)
    return "\t".join(
        cols[k]
        for k in (
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
        )
    )


def _paf_line(**kw) -> str:
    cols = {
        "qname": "q1",
        "qlen": "1000",
        "qstart": "100",
        "qend": "600",
        "strand": "+",
        "tname": "ref1",
        "tlen": "2000",
        "tstart": "300",
        "tend": "800",
        "nmatch": "490",
        "aln": "500",
        "mapq": "60",
    }
    cols.update(kw)
    base = "\t".join(
        cols[k]
        for k in (
            "qname",
            "qlen",
            "qstart",
            "qend",
            "strand",
            "tname",
            "tlen",
            "tstart",
            "tend",
            "nmatch",
            "aln",
            "mapq",
        )
    )
    return base + "\tNM:i:10"


class TestHitBlast:
    def test_full_line(self):
        h = Hit.from_blast_line(_blast_line())
        assert h.query_id == "q1"
        assert h.identity == 99.5
        assert h.query_coverage == pytest.approx(400 / 500 * 100)
        assert h.subject_coverage == pytest.approx(400 / 450 * 100)
        assert h.strand == "+"
        assert h.backend == "blast"
        assert h.mismatches == 2

    def test_minus_strand(self):
        h = Hit.from_blast_line(_blast_line(sstart="400", send="1"))
        assert h.strand == "-"

    def test_short_line_raises(self):
        with pytest.raises(ValueError):
            Hit.from_blast_line("a\tb\tc")

    def test_zero_lengths_safe(self):
        h = Hit.from_blast_line(_blast_line(qlen="0", slen="0"))
        assert h.query_coverage == 0.0 and h.subject_coverage == 0.0

    def test_lengths_populated(self):
        h = Hit.from_blast_line(_blast_line())
        assert h.query_length == 500 and h.subject_length == 450

    def test_to_dict_roundtrip(self):
        h = Hit.from_blast_line(_blast_line())
        assert h.to_dict()["query_id"] == "q1"


class TestHitPaf:
    def test_full_line(self):
        h = Hit.from_paf_line(_paf_line())
        assert h.identity == pytest.approx(490 / 500 * 100)
        assert h.mapq == 60
        assert h.mismatches == 10
        assert h.backend == "minimap2"

    def test_bad_nm_suppressed(self):
        h = Hit.from_paf_line(_paf_line().replace("NM:i:10", "NM:i:xx"))
        assert h.mismatches == 0

    def test_empty_comment_short_raise(self):
        for bad in ["", "# comment", "a\tb\tc"]:
            with pytest.raises(ValueError):
                Hit.from_paf_line(bad)

    def test_zero_aln_safe(self):
        h = Hit.from_paf_line(_paf_line(nmatch="0", aln="0"))
        assert h.identity == 0.0


def _hit(sstart: int, send: int, identity: float = 99.0, aln: int = 100) -> Hit:
    return Hit(subject_start=sstart, subject_end=send, identity=identity, alignment_length=aln)


class TestMergeIntervals:
    def test_empty(self):
        assert merge_intervals([]) == (0.0, 0.0)

    def test_overlap_merges(self):
        cov, avg = merge_intervals([_hit(1, 100), _hit(50, 150)], subject_length=200)
        assert cov == 75.0
        assert avg == 99.0

    def test_weighted_identity(self):
        cov, avg = merge_intervals(
            [_hit(1, 100, 100.0, 100), _hit(200, 300, 80.0, 100)],
            subject_length=400,
        )
        assert avg == 90.0

    def test_length_inferred(self):
        cov, _ = merge_intervals([_hit(1, 100)])
        assert cov == 100.0

    def test_zero_coords(self):
        assert merge_intervals([_hit(0, 0)]) == (0.0, 0.0)

    def test_cap_100(self):
        cov, _ = merge_intervals([_hit(1, 300)], subject_length=200)
        assert cov == 100.0

    def test_flipped_coords(self):
        cov, _ = merge_intervals([_hit(150, 50)], subject_length=200)
        assert cov == pytest.approx(101 / 200 * 100)


class TestConfidenceTier:
    @pytest.mark.parametrize(
        "cov,ident,expected",
        [
            (100.0, 100.0, "perfect"),
            (99.0, 99.0, "perfect"),
            (99.0, 95.0, "very_high"),
            (99.0, 94.9, "high"),
            (95.0, 90.0, "high"),
            (90.0, 85.0, "good"),
            (80.0, 80.0, "low"),
            (79.9, 100.0, "none"),
            (100.0, 79.9, "none"),
        ],
    )
    def test_boundaries(self, cov, ident, expected):
        assert confidence_tier(cov, ident) == expected


class TestClassifyAllele:
    def test_four_labels(self):
        assert classify_allele(100.0, 100.0) == ("exact", 90)
        assert classify_allele(96.0, 99.0) == ("novel", 63)
        assert classify_allele(96.0, 60.0) == ("partial", 18)
        assert classify_allele(70.0, 70.0) == ("missing", 0)

    def test_novel_boundary(self):
        assert classify_allele(95.0, 98.0)[0] == "novel"
        assert classify_allele(95.0, 97.9)[0] == "partial"


class TestMergeEdge:
    def test_no_valid_intervals_with_length(self):
        assert merge_intervals([_hit(0, 0, identity=99.0, aln=10)], subject_length=200) == (
            0.0,
            0.0,
        )
