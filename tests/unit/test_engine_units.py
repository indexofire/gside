"""Pure unit tests for engine/hits.py (no binaries)."""

from __future__ import annotations

import pytest

from gside.engine.hits import Hit


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
