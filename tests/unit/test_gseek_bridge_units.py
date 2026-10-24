"""Unit tests for the optional gseek bridge (no gseek import required)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from gside.engine import gseek_bridge as bridge


def test_engine_mode_defaults_to_legacy(monkeypatch):
    monkeypatch.delenv("GSIDE_ENGINE", raising=False)
    assert bridge.engine_mode() == "legacy"


@pytest.mark.parametrize("mode", ["legacy", "gseek", "shadow", " GSEEK "])
def test_engine_mode_accepts_known_values(monkeypatch, mode):
    monkeypatch.setenv("GSIDE_ENGINE", mode)
    assert bridge.engine_mode() == mode.strip().lower()


def test_engine_mode_rejects_unknown_value(monkeypatch):
    monkeypatch.setenv("GSIDE_ENGINE", "fastest")
    with pytest.raises(RuntimeError, match="unknown GSIDE_ENGINE"):
        bridge.engine_mode()


def _record(**overrides):
    values = {
        "query_id": "q1",
        "target_id": "markers~~~inva~~~x",
        "identity": 99.5,
        "aligned_length": 90,
        "query_length": 100,
        "target_length": 300,
        "query_start": 9,
        "query_end": 99,
        "target_start": 19,
        "target_end": 109,
        "evalue": 1e-5,
        "bit_score": 100.0,
        "is_reverse": False,
        "mismatches": 1,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_record_to_hit_preserves_gside_formulas_and_coordinates():
    hit = bridge.record_to_hit(_record())
    assert hit.query_id == "q1"
    assert hit.subject_id == "markers~~~inva~~~x"
    assert hit.identity == 99.5
    assert hit.query_coverage == pytest.approx(90.0)
    assert hit.subject_coverage == pytest.approx(30.0)
    assert (hit.query_start, hit.query_end) == (10, 99)
    assert (hit.subject_start, hit.subject_end) == (20, 109)
    assert hit.strand == "+"
    assert hit.backend == "blast"


def test_record_to_hit_maps_reverse_strand():
    hit = bridge.record_to_hit(_record(is_reverse=True))
    assert hit.strand == "-"


def test_aggregate_screen_unions_spans_and_weights_identity():
    records = [
        _record(target_start=0, target_end=150, aligned_length=150, matches=145),
        _record(target_start=100, target_end=250, aligned_length=150, matches=150),
    ]
    (hit,) = bridge.aggregate_screen(records)
    assert hit["gene"] == "inva"
    assert hit["coverage"] == pytest.approx(83.3)
    assert hit["identity"] == pytest.approx(98.33)
    assert hit["ref_len"] == 300
    assert hit["seqid"] == "markers~~~inva~~~x"


def test_aggregate_screen_skips_records_without_coordinates():
    records = [_record(target_start=None, target_end=None, aligned_length=150, matches=150)]
    assert bridge.aggregate_screen(records) == []


def test_aggregate_screen_uses_gene_fallback_for_plain_ids():
    records = [
        _record(
            target_id="plain-id",
            target_start=0,
            target_end=30,
            aligned_length=30,
            matches=30,
        )
    ]
    assert bridge.aggregate_screen(records)[0]["gene"] == "plain-id"


class _FakeGseek:
    def __init__(self, records):
        self.records = records
        self.calls = []

    def align(self, **kwargs):
        self.calls.append(kwargs)
        return SimpleNamespace(records=self.records)


def test_blast_hits_maps_gside_defaults_and_filters(monkeypatch):
    fake = _FakeGseek([_record(identity=99.5, aligned_length=90, query_length=100)])
    monkeypatch.setattr(bridge, "_require_gseek", lambda: fake)
    monkeypatch.setattr(bridge, "_resolver", lambda: "resolver")

    hits = bridge.blast_hits(
        "contigs.fna", "markers.fasta", evalue=1e-10, threads=2, min_identity=99.0
    )
    assert len(hits) == 1
    call = fake.calls[0]
    assert call["backend"] == "blast"
    assert call["target"] == "markers.fasta"
    assert call["threads"] == 2
    assert call["resolver"] == "resolver"
    assert call["params"] == {
        "task": "megablast",
        "evalue": 1e-10,
        "max_targets": 500,
        "word_size": 11,
        "dust": True,
    }


def test_blast_hits_filters_low_coverage(monkeypatch):
    fake = _FakeGseek([_record(identity=99.5, aligned_length=10, query_length=100)])
    monkeypatch.setattr(bridge, "_require_gseek", lambda: fake)
    monkeypatch.setattr(bridge, "_resolver", lambda: "resolver")
    assert bridge.blast_hits("q.fna", "m.fasta", min_coverage=50.0) == []


def test_screen_reads_gseek_passes_paired_fastq(monkeypatch):
    records = [_record(target_start=0, target_end=150, aligned_length=150, matches=145)]
    fake = _FakeGseek(records)
    monkeypatch.setattr(bridge, "_require_gseek", lambda: fake)
    monkeypatch.setattr(bridge, "_resolver", lambda: "resolver")

    hits = bridge.screen_reads_gseek(
        ["r1.fq", "r2.fq"], "markers.fasta", threads=8, preset="map-ont"
    )
    call = fake.calls[0]
    assert call["backend"] == "minimap2"
    assert call["query"] == (Path("r1.fq"), Path("r2.fq"))
    assert call["paired"] is True
    assert call["query_kind"] == "fastq"
    assert call["params"] == {"preset": "map-ont", "secondary": "no"}
    assert call["threads"] == 8
    assert hits[0]["gene"] == "inva"


def test_require_gseek_without_package_gives_actionable_error(monkeypatch):
    def missing(name):
        raise ImportError(name)

    monkeypatch.setattr(bridge, "importlib", SimpleNamespace(import_module=missing))
    with pytest.raises(RuntimeError, match="requires the optional gseek"):
        bridge._require_gseek()


def test_resolver_reuses_pixi_paths(monkeypatch, tmp_path):
    first = tmp_path / "bin1"
    second = tmp_path / "bin2"
    first.mkdir()
    second.mkdir()
    captured = {}

    class FakeResolver:
        def __init__(self, *, extra_paths):
            captured["extra_paths"] = extra_paths

    import gside.config as gside_config

    fake_execution = SimpleNamespace(ToolResolver=FakeResolver)
    monkeypatch.setattr(
        bridge, "importlib", SimpleNamespace(import_module=lambda name: fake_execution)
    )
    monkeypatch.setattr(gside_config, "pixi_path", lambda: f"{first}:{second}:/does/not/exist")
    bridge._resolver()
    assert captured["extra_paths"] == [first, second]


def test_gseek_available_uses_find_spec(monkeypatch):
    monkeypatch.setattr(bridge.importlib.util, "find_spec", lambda name: object())
    assert bridge.gseek_available() is True


class _FakeDistanceGseek:
    def __init__(self, records):
        self.records = records
        self.calls = []

    def distance(self, **kwargs):
        self.calls.append(kwargs)
        return SimpleNamespace(records=self.records)


def test_mash_distance_maps_external_reference(monkeypatch):
    record = SimpleNamespace(
        query_id="query",
        reference_id="GCF_1",
        distance=0.01,
        pvalue=0.0,
        shared_hashes=990,
        total_hashes=1000,
    )
    fake = _FakeDistanceGseek([record])
    monkeypatch.setattr(bridge, "_require_gseek", lambda: fake)
    monkeypatch.setattr(bridge, "_resolver", lambda: "resolver")

    hits = bridge.mash_distance("q.fna", "refs.msh", threads=2, max_distance=1.0)
    call = fake.calls[0]
    assert call["backend"] == "mash"
    assert call["reference"] == Path("refs.msh")
    assert call["external_reference"] is True
    assert call["params"] == {"max_distance": 1.0}
    assert hits[0].reference_id == "GCF_1"
    assert hits[0].distance == 0.01


def test_sourmash_gather_maps_rows(monkeypatch):
    record = SimpleNamespace(
        reference_id="abc123",
        ani=99.0,
        tags={"reference_name": "GCF_1 Escherichia coli", "f_unique_weighted": "0.42"},
    )
    fake = _FakeDistanceGseek([record])
    monkeypatch.setattr(bridge, "_require_gseek", lambda: fake)
    monkeypatch.setattr(bridge, "_resolver", lambda: "resolver")

    rows = bridge.sourmash_gather("q.fna", "gtdb.zip")
    call = fake.calls[0]
    assert call["backend"] == "sourmash"
    assert call["external_reference"] is True
    assert call["params"]["mode"] == "gather"
    assert rows == [
        {
            "name": "GCF_1 Escherichia coli",
            "f_unique_weighted": "0.42",
            "average_containment_ani": "0.99",
        }
    ]


def test_sourmash_gather_rejects_multiple_query_files(monkeypatch):
    monkeypatch.setattr(bridge, "_require_gseek", lambda: object())
    with pytest.raises(RuntimeError, match="one query file"):
        bridge.sourmash_gather(["a.fna", "b.fna"], "gtdb.zip")


def test_log_mash_shadow_accepts_legacy_tuples_and_records(caplog):
    from gside.engine.backends.kmer import KmerDistance

    with caplog.at_level("WARNING"):
        bridge.log_mash_shadow(
            [("GCF_1", 0.99)],
            [KmerDistance("q", "GCF_2", 0.01, 0.0, 990, 1000)],
        )
    assert "gseek shadow mash differs" in caplog.text
