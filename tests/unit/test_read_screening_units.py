"""Unit tests for minimap2 read screening (PAF parsing + breadth aggregation)."""

from __future__ import annotations

import gzip

import pytest

from gside.analysis import read_screening as rs


class TestIsFastq:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("r.fq", True),
            ("r.fastq", True),
            ("r.fq.gz", True),
            ("r.FASTQ.GZ", True),
            ("g.fna", False),
            ("g.fasta", False),
            ("notes.txt", False),
        ],
    )
    def test_suffixes(self, name, expected):
        assert rs.is_fastq(name) is expected


class TestMeanReadLength:
    def test_plain_fastq(self, tmp_path):
        fq = tmp_path / "r.fq"
        fq.write_text("@a\nACGT\n+\nIIII\n@b\nACGTACGT\n+\nIIIIIIII\n")
        assert rs.mean_read_length(fq) == pytest.approx(6.0)

    def test_gzip_fastq(self, tmp_path):
        fq = tmp_path / "r.fq.gz"
        with gzip.open(fq, "wt") as fh:
            fh.write("@a\nACGTACGTAC\n+\nIIIIIIIIII\n")
        assert rs.mean_read_length(fq) == pytest.approx(10.0)

    def test_empty_is_zero(self, tmp_path):
        fq = tmp_path / "r.fq"
        fq.write_text("")
        assert rs.mean_read_length(fq) == 0.0


class TestPafParsing:
    def test_primary_row(self):
        row = rs._parse_paf_line(
            "r1\t150\t0\t150\t+\tmarkers~~~inva~~~M90846.1\t900\t10\t160\t148\t150\t60\ttp:A:P"
        )
        assert row == {
            "qname": "r1",
            "tname": "markers~~~inva~~~M90846.1",
            "tlen": 900,
            "tstart": 10,
            "tend": 160,
            "nmatch": 148,
            "alen": 150,
        }

    def test_secondary_row_skipped(self):
        assert (
            rs._parse_paf_line(
                "r1\t150\t0\t150\t+\tmarkers~~~inva~~~x\t900\t10\t160\t148\t150\t0\ttp:A:S"
            )
            is None
        )

    def test_short_row_skipped(self):
        assert rs._parse_paf_line("r1\t150\t0") is None

    def test_non_numeric_row_skipped(self):
        assert rs._parse_paf_line("a\tb\tc\td\te\tf\tg\th\ti\tj\tk") is None


class TestAggregate:
    def _row(self, tstart, tend, nmatch=145, alen=150, tlen=300, tname="markers~~~inva~~~x"):
        return {
            "qname": "r",
            "tname": tname,
            "tlen": tlen,
            "tstart": tstart,
            "tend": tend,
            "nmatch": nmatch,
            "alen": alen,
        }

    def test_breadth_is_union_of_intervals(self):
        hits = rs._aggregate([self._row(0, 150), self._row(100, 250)])
        assert hits[0]["coverage"] == pytest.approx(83.3)
        assert hits[0]["ref_len"] == 300
        assert hits[0]["gene"] == "inva"

    def test_identity_is_alignment_weighted(self):
        hits = rs._aggregate(
            [
                self._row(0, 100, nmatch=100, alen=100),
                self._row(100, 200, nmatch=80, alen=100),
            ]
        )
        assert hits[0]["identity"] == pytest.approx(90.0)

    def test_multiple_targets(self):
        hits = rs._aggregate(
            [
                self._row(0, 150, tname="markers~~~inva~~~x"),
                self._row(0, 60, tname="markers~~~mapa~~~y", tlen=120),
            ]
        )
        by_gene = {hit["gene"]: hit for hit in hits}
        assert by_gene["inva"]["coverage"] == pytest.approx(50.0)
        assert by_gene["mapa"]["coverage"] == pytest.approx(50.0)

    def test_zero_length_alignment_skipped(self):
        assert rs._aggregate([self._row(0, 0, alen=0)]) == []


class TestScreenReads:
    def _stub(self, monkeypatch, stdout, mean=150.0):
        seen = {}

        def fake_run(cmd, **kwargs):
            seen["cmd"] = cmd
            return type("P", (), {"stdout": stdout})()

        monkeypatch.setattr(rs, "require_bin", lambda *a, **k: "/bin/minimap2")
        monkeypatch.setattr(rs, "run_checked", fake_run)
        monkeypatch.setattr(rs, "mean_read_length", lambda path: mean)
        return seen

    def test_short_reads_use_sr_preset(self, monkeypatch, tmp_path):
        seen = self._stub(monkeypatch, "")
        rs.screen_reads("r.fq", "markers.fasta", threads=8)
        cmd = seen["cmd"]
        assert cmd[1] == "-x" and cmd[2] == "sr"
        assert "--secondary=no" in cmd
        assert cmd[cmd.index("-t") + 1] == "8"
        assert cmd[-2:] == ["markers.fasta", "r.fq"]

    def test_long_reads_use_map_ont(self, monkeypatch, tmp_path):
        seen = self._stub(monkeypatch, "", mean=5000.0)
        rs.screen_reads("r.fq", "markers.fasta")
        assert seen["cmd"][2] == "map-ont"

    def test_rows_are_aggregated(self, monkeypatch):
        paf = (
            "r1\t150\t0\t150\t+\tmarkers~~~inva~~~x\t300\t0\t150\t145\t150\t60\ttp:A:P\n"
            "r2\t150\t0\t150\t+\tmarkers~~~inva~~~x\t300\t150\t300\t145\t150\t60\ttp:A:P\n"
        )
        self._stub(monkeypatch, paf)
        hits = rs.screen_reads("r.fq", "markers.fasta")
        assert hits[0]["gene"] == "inva"
        assert hits[0]["coverage"] == pytest.approx(100.0)

    def test_no_alignments_is_empty(self, monkeypatch):
        self._stub(monkeypatch, "")
        assert rs.screen_reads("r.fq", "markers.fasta") == []


class TestEngineSelection:
    def test_gseek_delegates_to_bridge(self, monkeypatch):
        monkeypatch.setenv("GSIDE_ENGINE", "gseek")
        monkeypatch.setattr(rs, "mean_read_length", lambda path: 150.0)
        seen = {}

        def fake_screen(files, markers, *, threads, preset):
            seen.update(files=files, markers=markers, threads=threads, preset=preset)
            return [{"gene": "inva"}]

        monkeypatch.setattr(rs, "screen_reads_gseek", fake_screen)
        assert rs.screen_reads("r.fq", "markers.fasta", threads=2) == [{"gene": "inva"}]
        assert seen == {
            "files": ["r.fq"],
            "markers": "markers.fasta",
            "threads": 2,
            "preset": "sr",
        }

    def test_shadow_runs_both_and_returns_legacy(self, monkeypatch):
        monkeypatch.setenv("GSIDE_ENGINE", "shadow")
        monkeypatch.setattr(rs, "mean_read_length", lambda path: 150.0)
        monkeypatch.setattr(
            rs, "_screen_reads_legacy", lambda files, markers, threads, preset: [{"gene": "old"}]
        )
        monkeypatch.setattr(
            rs,
            "screen_reads_gseek",
            lambda files, markers, *, threads, preset: [{"gene": "new"}],
        )
        seen = {}
        monkeypatch.setattr(rs, "log_screen_shadow", lambda old, new: seen.update(old=old, new=new))
        assert rs.screen_reads("r.fq", "markers.fasta") == [{"gene": "old"}]
        assert seen == {"old": [{"gene": "old"}], "new": [{"gene": "new"}]}

    def test_invalid_engine_mode_fails_before_tool_discovery(self, monkeypatch):
        monkeypatch.setenv("GSIDE_ENGINE", "fastest")
        monkeypatch.setattr(rs, "mean_read_length", lambda path: 150.0)
        with pytest.raises(RuntimeError, match="unknown GSIDE_ENGINE"):
            rs.screen_reads("r.fq", "markers.fasta")


class TestMalformedInputs:
    def test_non_at_header_is_skipped(self, tmp_path):
        fq = tmp_path / "r.fq"
        fq.write_text("junk line\n@a\nACGT\n+\nIIII\n")
        assert rs.mean_read_length(fq) == pytest.approx(4.0)

    def test_paf_non_numeric_field_skipped(self):
        row = "r1\t150\t0\t150\t+\tt\tnotanint\t0\t150\t145\t150\t60\ttp:A:P"
        assert rs._parse_paf_line(row) is None

    def test_contained_interval_is_ignored(self):
        row = {
            "qname": "r",
            "tname": "markers~~~x~~~y",
            "tlen": 300,
            "tstart": 0,
            "tend": 300,
            "nmatch": 300,
            "alen": 300,
        }
        inner = dict(row, tstart=10, tend=50, nmatch=40, alen=40)
        hits = rs._aggregate([row, inner])
        assert hits[0]["coverage"] == pytest.approx(100.0)
