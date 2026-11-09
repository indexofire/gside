"""Multi-database runs for the ANI and sourmash identifiers (stubs only)."""

from __future__ import annotations

from pathlib import Path

import pytest

import gside.analysis.ani_identifier as ani
import gside.analysis.sourmash_identifier as sm
from gside.analysis._envelope import MethodEnvelope
from gside.analysis.ani_identifier import AniHit, identify_by_ani


def _skani_hit(ref="GCF_A", ani=99.0, af=0.9):
    return AniHit(ref=ref, ani=ani, aligned_fraction=af)


@pytest.fixture
def taxadb(tmp_path: Path) -> Path:
    (tmp_path / "metadata.tsv").write_text("genome\tspecies\nGCF_A\tEscherichia coli K12\n")
    return tmp_path


@pytest.fixture
def salmdb(tmp_path: Path) -> Path:
    db = tmp_path / "lab"
    db.mkdir()
    (db / "metadata.tsv").write_text("genome\tspecies\nGCF_B\tSalmonella enterica\n")
    return db


class TestPanelMultiDatabase:
    def test_best_hit_across_databases(self, monkeypatch, taxadb, salmdb):
        hits_by_db = {
            str(taxadb): [_skani_hit(ref="GCF_A", ani=94.0)],
            str(salmdb): [_skani_hit(ref="GCF_B", ani=99.5)],
        }
        monkeypatch.setattr(ani, "_skani_search", lambda q, d, threads=4: hits_by_db[str(d)])
        identified = identify_by_ani("q.fna", mode="ani", dbs=[taxadb, salmdb])
        r = identified.result
        assert (r["species"], r["confidence"]) == ("Salmonella", "high")
        assert r["top_hits"][0]["genome"] == "GCF_B"
        assert identified.database["name"] == f"{taxadb.name}+{salmdb.name}"
        assert [d["name"] for d in identified.to_dict()["databases"]] == [
            taxadb.name,
            salmdb.name,
        ]

    def test_taxa_maps_merge_first_wins(self, monkeypatch, tmp_path):
        first = tmp_path / "one"
        first.mkdir()
        (first / "metadata.tsv").write_text("genome\tspecies\nGCF_A\tEscherichia coli K12\n")
        second = tmp_path / "two"
        second.mkdir()
        (second / "metadata.tsv").write_text("genome\tspecies\nGCF_A\tSalmonella enterica\n")
        monkeypatch.setattr(ani, "_skani_search", lambda q, d, threads=4: [_skani_hit(ref="GCF_A")])
        r = identify_by_ani("q.fna", mode="ani", dbs=[first, second]).result
        assert r["species"] == "Escherichia_coli"

    def test_single_db_list_keeps_single_shape(self, monkeypatch, taxadb):
        monkeypatch.setattr(ani, "_skani_search", lambda q, d, threads=4: [_skani_hit()])
        identified = identify_by_ani("q.fna", mode="ani", dbs=[taxadb])
        assert identified.database["name"] == taxadb.name
        assert "databases" not in identified.to_dict()


class TestMashMultiDatabase:
    def test_best_identity_across_databases(self, monkeypatch, tmp_path):
        one = tmp_path / "one"
        one.mkdir()
        (one / "mash.msh").write_bytes(b"msh")
        two = tmp_path / "two"
        two.mkdir()
        (two / "mash.msh").write_bytes(b"msh")
        rows_by_msh = {
            str(one / "mash.msh"): [("GCF_A", 0.95)],
            str(two / "mash.msh"): [("GCF_B", 0.995)],
        }
        monkeypatch.setattr(ani, "_mash_dist", lambda q, msh, threads=4: rows_by_msh[str(msh)])
        identified = identify_by_ani("q.fna", mode="mash", dbs=[one, two])
        r = identified.result
        assert r["identity"] == pytest.approx(0.995)
        assert r["top_hits"][0]["genome"] == "GCF_B"
        assert identified.database["name"] == "one+two"

    def test_reads_use_dist_reads_per_db(self, monkeypatch, tmp_path):
        one = tmp_path / "one"
        one.mkdir()
        (one / "mash.msh").write_bytes(b"msh")
        two = tmp_path / "two"
        two.mkdir()
        (two / "mash.msh").write_bytes(b"msh")
        rows_by_msh = {
            str(one / "mash.msh"): [("GCF_A", 0.9)],
            str(two / "mash.msh"): [("GCF_B", 0.98)],
        }
        monkeypatch.setattr(
            ani, "_mash_dist_reads", lambda files, msh, threads=4: rows_by_msh[str(msh)]
        )
        reads = tmp_path / "r.fq"
        reads.write_text("@r\nACGT\n+\nI\n")
        r = identify_by_ani(str(reads), mode="mash", dbs=[one, two]).result
        assert r["top_hits"][0]["genome"] == "GCF_B"


class TestSourmashMultiDatabase:
    def test_partition_merged_and_reranked(self, monkeypatch, tmp_path):
        one = tmp_path / "one"
        one.mkdir()
        (one / "gtdb-reps-k31.zip").write_bytes(b"zip")
        two = tmp_path / "two"
        two.mkdir()
        (two / "gtdb-reps-k31.zip").write_bytes(b"zip")
        partitions = {
            one: [
                {"species": "Escherichia coli", "f_unique_weighted": 0.5, "ani": 96.0},
            ],
            two: [
                {"species": "Salmonella enterica", "f_unique_weighted": 0.7, "ani": 96.5},
            ],
        }
        monkeypatch.setattr(sm, "_gather_and_tax", lambda contigs, base: partitions[base])
        identified = sm.identify_by_sourmash("q.fna", dbs=[one, two])
        r = identified.result
        assert r["f_unique_weighted"] == pytest.approx(0.7)
        assert identified.database["name"] == "one+two"
        assert [d["name"] for d in identified.to_dict()["databases"]] == ["one", "two"]

    def test_single_db_list_keeps_single_shape(self, monkeypatch, tmp_path):
        one = tmp_path / "one"
        one.mkdir()
        (one / "gtdb-reps-k31.zip").write_bytes(b"zip")
        monkeypatch.setattr(
            sm,
            "_gather_and_tax",
            lambda contigs, base: [
                {"species": "Escherichia coli", "f_unique_weighted": 0.5, "ani": 96.0}
            ],
        )
        identified = sm.identify_by_sourmash("q.fna", dbs=[one])
        assert identified.database["name"] == "one"
        assert "databases" not in identified.to_dict()


def test_envelope_databases_field_conditional():
    base = MethodEnvelope(method="panel", database={"name": "a"}, result={})
    assert "databases" not in base.to_dict()
    multi = MethodEnvelope(
        method="panel",
        database={"name": "a+b"},
        result={},
        databases=[{"name": "a"}, {"name": "b"}],
    )
    assert multi.to_dict()["databases"] == [{"name": "a"}, {"name": "b"}]


class TestMashMultiBands:
    def test_medium_band_across_databases(self, monkeypatch, tmp_path):
        one = tmp_path / "one"
        one.mkdir()
        (one / "mash.msh").write_bytes(b"msh")
        two = tmp_path / "two"
        two.mkdir()
        (two / "mash.msh").write_bytes(b"msh")
        rows_by_msh = {
            str(one / "mash.msh"): [("GCF_A", 0.9)],
            str(two / "mash.msh"): [("GCF_B", 0.94)],
        }
        monkeypatch.setattr(ani, "_mash_dist", lambda q, msh, threads=4: rows_by_msh[str(msh)])
        r = identify_by_ani("q.fna", mode="mash", dbs=[one, two]).result
        assert r["confidence"] == "medium"
