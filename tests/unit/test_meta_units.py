"""Pure unit tests for gside meta: kraken2/sylph identifiers + `gside meta` CLI (no binaries)."""

from __future__ import annotations

import json
import types
from pathlib import Path

import pytest
from click.testing import CliRunner

import gside.analysis.kraken2_identifier as k2
import gside.analysis.sylph_identifier as sy
from gside import __version__, cli
from gside.cli import _run_meta, main

# ---------------------------------------------------------------- kraken2 fakes

_KREPORT = "\n".join(
    [
        "  1.00\t100\t100\tR\t1\troot",
        "  0.80\t80\t80\tS\t562\t  Escherichia coli",
        " 99.00\t9900\t9900\tU\t0\tunclassified",
    ]
)

_BRACKEN = (
    "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads"
    "\tfraction_total_reads\n"
    "Salmonella enterica\t289\tS\t40\t10\t50\t0.0200\n"
    "Escherichia coli\t562\tS\t800\t50\t850\t0.8500\n"
    "Campylobacter jejuni\t197\tS\t1\t0\t1\t0.0010\n"
    "Vibrio sp.\t10\tspecies\t5\t0\t5\t0.0100\n"
    "Enterobacteriaceae\t543\tF\t10\t0\t10\t0.1000\n"
)


def _k2_run(
    kreport: str = _KREPORT,
    bracken: str = _BRACKEN,
    cmds: list | None = None,
    fail_tool: str | None = None,
):
    """Fake run_checked mimicking kraken2 (writes kreport) then bracken (writes bracken.txt)."""

    def run_checked(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
        if cmds is not None:
            cmds.append((list(cmd), timeout, name))
        if fail_tool and cmd[0].endswith(fail_tool):
            raise RuntimeError(f"{name} (exit 1): boom")
        work = (
            Path(cmd[cmd.index("--report") + 1]).parent
            if "--report" in cmd
            else Path(cmd[cmd.index("-i") + 1]).parent
        )
        if "--report" in cmd:  # kraken2 invocation
            (work / "kraken.out").write_text("")
            (work / "kreport.txt").write_text(kreport)
        else:  # bracken invocation
            (work / "bracken.txt").write_text(bracken)
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    return run_checked


def _k2_ready(monkeypatch, tmp_path: Path) -> Path:
    db = tmp_path / "D5_kraken2"
    db.mkdir()
    (db / "hash.k2d").write_bytes(b"x")
    monkeypatch.setattr(k2, "which", lambda name: f"/bin/{name}")
    return db


# ---------------------------------------------------------------- sylph fakes

_PROFILE = (
    "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov\n"
    "sample.fq\tEscherichia_coli.fna\t78.12\t81.82\t97.53\t264.000\n"
    "sample.fq\tSalmonella_enterica.fna\t10.50\t11.00\t94.10\t30.200\n"
    "sample.fq\tStaphylococcus_aureus.fna\t5.00\t4.00\t91.00\t12.000\n"
)


def _sy_run(profile: str = _PROFILE, cmds: list | None = None, sketch_suffix: str = ".sylsp"):
    """Fake run_checked mimicking sylph sketch then profile (writes profile.tsv)."""

    def run_checked(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
        if cmds is not None:
            cmds.append((list(cmd), timeout, name))
        if cmd[1] == "sketch":
            # v1 writes <sketches dir>/<sample>.sylsp (or .syldata for v0)
            sketch_dir = Path(cmd[cmd.index("-d") + 1])
            sketch_dir.mkdir(parents=True, exist_ok=True)
            (sketch_dir / f"reads.fq{sketch_suffix}").write_bytes(b"x")
        else:  # profile invocation: output path follows -o
            Path(cmd[cmd.index("-o") + 1]).write_text(profile)
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    return run_checked


def _sy_ready(monkeypatch, tmp_path: Path) -> Path:
    db = tmp_path / "D6_sylph"
    db.mkdir()
    (db / "gtdb-r226.syldb").write_bytes(b"x")
    monkeypatch.setattr(sy, "which", lambda name: "/bin/sylph")
    return db


# ================================================================ kraken2 module


class TestKraken2ToDict:
    def test_shape(self):
        d = k2.Kraken2Result(database={"name": "D5_kraken2"}, result={"top_species": []}).to_dict()
        assert d["analysis_type"] == "metagenomic_profiling"
        assert d["method"] == "kraken2_bracken"
        assert d["database"] == {"name": "D5_kraken2"}
        assert d["result"] == {"top_species": []}


class TestKraken2Unavailable:
    def test_missing_binary_raises_actionable(self, monkeypatch, tmp_path):
        db = tmp_path / "D5_kraken2"
        db.mkdir()
        (db / "hash.k2d").write_bytes(b"x")
        monkeypatch.setattr(k2, "which", lambda name: None)
        with pytest.raises(RuntimeError, match="kraken2 not found") as exc:
            k2.identify_by_kraken2(["r.fq"], db_dir=str(db))
        assert "pixi add kraken2 bracken" in str(exc.value)

    def test_missing_db_raises_actionable(self, monkeypatch, tmp_path):
        monkeypatch.setattr(k2, "which", lambda name: f"/bin/{name}")
        with pytest.raises(RuntimeError, match="Kraken2 database not found") as exc:
            k2.identify_by_kraken2(["r.fq"], db_dir=str(tmp_path / "empty"))
        assert "gside db setup --tier meta" in str(exc.value)

    def test_missing_db_default_path(self, monkeypatch, tmp_path):
        monkeypatch.setattr(k2, "which", lambda name: f"/bin/{name}")
        monkeypatch.setattr(k2, "SPECIES_DB_DIR", tmp_path)
        with pytest.raises(RuntimeError, match="Kraken2 database not found"):
            k2.identify_by_kraken2(["r.fq"])


class TestKraken2HappyPath:
    def test_species_parsed_sorted_with_confidence(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        monkeypatch.setattr(k2, "run_checked", _k2_run())
        r = k2.identify_by_kraken2(["r.fq"], db_dir=str(db))
        species = r.result["top_species"]
        assert [s["name"] for s in species] == [
            "Escherichia coli",
            "Salmonella enterica",
            "Vibrio sp.",
            "Campylobacter jejuni",
        ]
        assert species[0] == {
            "name": "Escherichia coli",
            "taxid": 562,
            "reads": 850,
            "fraction": 0.85,
            "confidence": "high",
        }
        assert [s["confidence"] for s in species] == ["high", "medium", "medium", "low"]

    def test_non_species_rows_filtered(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        monkeypatch.setattr(k2, "run_checked", _k2_run())
        r = k2.identify_by_kraken2(["r.fq"], db_dir=str(db))
        assert all(s["name"] != "Enterobacteriaceae" for s in r.result["top_species"])
        assert r.result["n_species"] == 4

    def test_totals_from_kreport(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        monkeypatch.setattr(k2, "run_checked", _k2_run())
        r = k2.identify_by_kraken2(["r.fq"], db_dir=str(db))
        assert r.result["total_reads"] == 10000
        assert r.result["classified_reads"] == 100

    def test_top_species_capped_at_20(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        rows = "".join(f"Species {i}\t{i}\tS\t1\t0\t1\t{0.5 - i * 0.01:.4f}\n" for i in range(25))
        header = _BRACKEN.splitlines()[0] + "\n"
        monkeypatch.setattr(k2, "run_checked", _k2_run(bracken=header + rows))
        r = k2.identify_by_kraken2(["r.fq"], db_dir=str(db))
        assert len(r.result["top_species"]) == 20
        assert r.result["n_species"] == 25

    def test_command_shape_paired(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        cmds: list = []
        monkeypatch.setattr(k2, "run_checked", _k2_run(cmds=cmds))
        k2.identify_by_kraken2(["r1.fq", "r2.fq"], db_dir=str(db), threads=9)
        k2_cmd, k2_to, k2_name = cmds[0]
        assert k2_cmd[0] == "/bin/kraken2"
        assert k2_cmd[k2_cmd.index("--db") + 1] == str(db)
        assert k2_cmd[k2_cmd.index("--threads") + 1] == "9"
        assert k2_cmd[-2:] == ["r1.fq", "r2.fq"]
        assert "--paired" in k2_cmd
        assert k2_to == 3600 and k2_name == "kraken2 failed"
        br_cmd, br_to, br_name = cmds[1]
        assert br_cmd[0] == "/bin/bracken"
        assert br_cmd[br_cmd.index("-d") + 1] == str(db)
        assert br_cmd[br_cmd.index("-r") + 1] == "150"
        assert br_cmd[br_cmd.index("-t") + 1] == "10"
        assert br_cmd[br_cmd.index("-o") + 1].endswith("bracken.txt")
        assert br_to == 300 and br_name == "bracken failed"

    def test_command_shape_single_has_no_paired_flag(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        cmds: list = []
        monkeypatch.setattr(k2, "run_checked", _k2_run(cmds=cmds))
        k2.identify_by_kraken2(["contigs.fna"], db_dir=str(db))
        assert "--paired" not in cmds[0][0]

    def test_malformed_bracken_rows_skipped(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        header = _BRACKEN.splitlines()[0] + "\n"
        bracken = header + "Short row\t1\tS\nBroken taxid\tabc\tS\t1\t0\t1\t0.5\n"
        monkeypatch.setattr(k2, "run_checked", _k2_run(bracken=bracken))
        r = k2.identify_by_kraken2(["r.fq"], db_dir=str(db))
        assert r.result["top_species"] == []
        assert r.result["n_species"] == 0

    def test_short_kreport_lines_skipped(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        kreport = "junk line\n  1.00\t100\t100\tR\t1\troot\n"
        monkeypatch.setattr(k2, "run_checked", _k2_run(kreport=kreport))
        r = k2.identify_by_kraken2(["r.fq"], db_dir=str(db))
        assert r.result["total_reads"] == 100
        assert r.result["classified_reads"] == 100

    def test_kraken2_failure_propagates(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        monkeypatch.setattr(k2, "run_checked", _k2_run(fail_tool="kraken2"))
        with pytest.raises(RuntimeError, match="kraken2 failed"):
            k2.identify_by_kraken2(["r.fq"], db_dir=str(db))

    def test_bracken_failure_propagates(self, monkeypatch, tmp_path):
        db = _k2_ready(monkeypatch, tmp_path)
        monkeypatch.setattr(k2, "run_checked", _k2_run(fail_tool="bracken"))
        with pytest.raises(RuntimeError, match="bracken failed"):
            k2.identify_by_kraken2(["r.fq"], db_dir=str(db))


# ================================================================ sylph module


class TestSylphToDict:
    def test_shape(self):
        d = sy.SylphResult(database={"name": "D6_sylph"}, result={"species": []}).to_dict()
        assert d["analysis_type"] == "metagenomic_profiling"
        assert d["method"] == "sylph_ani"
        assert d["result"] == {"species": []}


class TestSylphSpeciesName:
    def test_strips_fna_and_underscores(self):
        assert sy._species_from_genome_file("Escherichia_coli.fna") == "Escherichia coli"

    def test_path_only_keeps_filename(self):
        assert sy._species_from_genome_file("/db/genomes/Salmonella_enterica.fa") == (
            "Salmonella enterica"
        )


class TestSylphUnavailable:
    def test_missing_binary_raises_actionable(self, monkeypatch, tmp_path):
        db = tmp_path / "D6_sylph"
        db.mkdir()
        (db / "db.syldb").write_bytes(b"x")
        monkeypatch.setattr(sy, "which", lambda name: None)
        with pytest.raises(RuntimeError, match="sylph not found") as exc:
            sy.identify_by_sylph(["r.fq"], db_dir=str(db))
        assert "pixi add sylph" in str(exc.value)

    def test_fasta_input_rejected(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)
        with pytest.raises(RuntimeError, match="sylph requires FASTQ reads") as exc:
            sy.identify_by_sylph(["contigs.fna"], db_dir=str(db))
        assert "--mode kraken2" in str(exc.value)

    def test_missing_db_raises_actionable(self, monkeypatch, tmp_path):
        monkeypatch.setattr(sy, "which", lambda name: "/bin/sylph")
        with pytest.raises(RuntimeError, match="sylph database not found") as exc:
            sy.identify_by_sylph(["r.fq"], db_dir=str(tmp_path / "empty"))
        assert "gside db setup --tier meta_sylph" in str(exc.value)

    def test_missing_db_default_path(self, monkeypatch, tmp_path):
        monkeypatch.setattr(sy, "which", lambda name: "/bin/sylph")
        monkeypatch.setattr(sy, "SPECIES_DB_DIR", tmp_path)
        with pytest.raises(RuntimeError, match="sylph database not found"):
            sy.identify_by_sylph(["r.fq"])


class TestSylphNameResolution:
    def test_metadata_tsv_maps_full_name_and_stem(self, tmp_path):
        (tmp_path / "metadata.tsv").write_text(
            "genome_path\tspecies\n"
            "/genomes/GCF_000005845.2_ASM584v2_genomic.fna\tEscherichia coli\n"
        )
        mapping = sy._load_taxa_map(tmp_path)
        assert mapping["/genomes/GCF_000005845.2_ASM584v2_genomic.fna"] == "Escherichia coli"
        assert mapping["GCF_000005845.2_ASM584v2_genomic"] == "Escherichia coli"

    def test_missing_metadata_returns_empty(self, tmp_path):
        assert sy._load_taxa_map(tmp_path / "nested") == {}

    def test_short_rows_skipped(self, tmp_path):
        (tmp_path / "metadata.tsv").write_text("header\nonlyonecolumn\n")
        assert sy._load_taxa_map(tmp_path) == {}

    def test_genome_stem_fallback(self):
        assert sy._species_from_genome_file("/x/Escherichia_coli.fna") == "Escherichia coli"
        assert sy._species_from_genome_file("/x/GCF_000005845.2.fna") == "GCF_000005845.2"


class TestSylphHappyPath:
    def test_profile_parsed_with_ani_confidence(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)
        monkeypatch.setattr(sy, "run_checked", _sy_run())
        r = sy.identify_by_sylph(["r.fq"], db_dir=str(db))
        assert r.database == {"name": "D6_sylph", "path": str(db / "gtdb-r226.syldb")}
        assert r.result["n_species"] == 3
        first = r.result["species"][0]
        assert first["name"] == "Escherichia coli"
        assert first["ani"] == 97.53
        assert first["abundance"] == pytest.approx(0.7812)
        assert first["coverage"] == 264.0
        assert first["confidence"] == "high"
        assert [s["confidence"] for s in r.result["species"]] == ["high", "medium", "low"]

    def test_command_shape(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)
        cmds: list = []
        monkeypatch.setattr(sy, "run_checked", _sy_run(cmds=cmds))
        sy.identify_by_sylph(["r1.fq", "r2.fq"], db_dir=str(db), threads=7)
        sk_cmd, sk_to, sk_name = cmds[0]
        assert sk_cmd[:3] == ["/bin/sylph", "sketch", "-t"]
        assert sk_cmd[3] == "7"
        # paired files ride the -1/-2 flags (sylph v1); sketch dir via -d
        assert sk_cmd[sk_cmd.index("-1") + 1] == "r1.fq"
        assert sk_cmd[sk_cmd.index("-2") + 1] == "r2.fq"
        assert sk_cmd[sk_cmd.index("-d") + 1].endswith("sketches")
        assert sk_to == 3600 and sk_name == "sylph sketch failed"
        pr_cmd, pr_to, pr_name = cmds[1]
        assert pr_cmd[:3] == ["/bin/sylph", "profile", "-d"]
        assert pr_cmd[3] == str(db / "gtdb-r226.syldb")
        assert pr_cmd[4].endswith(".sylsp")
        assert pr_cmd[pr_cmd.index("-o") + 1].endswith("profile.tsv")
        assert pr_to == 600 and pr_name == "sylph profile failed"

    def test_single_end_sketch_passes_files_positionally(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)
        cmds: list = []
        monkeypatch.setattr(sy, "run_checked", _sy_run(cmds=cmds))
        sy.identify_by_sylph(["r.fq"], db_dir=str(db))
        sk_cmd = cmds[0][0]
        assert "-1" not in sk_cmd and "-2" not in sk_cmd
        assert "r.fq" in sk_cmd

    def test_sketch_output_discovered_across_suffixes(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)
        for suffix in (".sylsp", ".syldata"):
            cmds: list = []
            monkeypatch.setattr(sy, "run_checked", _sy_run(cmds=cmds, sketch_suffix=suffix))
            sy.identify_by_sylph(["r.fq"], db_dir=str(db))
            profile_cmd = cmds[1][0]
            sketch_arg = profile_cmd[4]
            assert Path(sketch_arg).suffix in {".sylsp", ".syldata"}

    def test_malformed_profile_rows_skipped(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)
        profile = (
            "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI"
            "\tEff_cov\n"
            "s.fq\tBad_values.fna\tNA\tNA\tNA\tNA\n"
        )
        monkeypatch.setattr(sy, "run_checked", _sy_run(profile=profile))
        r = sy.identify_by_sylph(["r.fq"], db_dir=str(db))
        assert r.result["species"] == []
        assert r.result["n_species"] == 0

    def test_header_only_profile_is_empty(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)
        header = _PROFILE.splitlines()[0] + "\n"
        monkeypatch.setattr(sy, "run_checked", _sy_run(profile=header))
        r = sy.identify_by_sylph(["r.fq"], db_dir=str(db))
        assert r.result["species"] == []

    def test_sketch_writing_nothing_raises(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)

        def silent_run(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(sy, "run_checked", silent_run)
        with pytest.raises(RuntimeError, match="sylph sketch produced no"):
            sy.identify_by_sylph(["r.fq"], db_dir=str(db))

    def test_sketch_failure_propagates(self, monkeypatch, tmp_path):
        db = _sy_ready(monkeypatch, tmp_path)

        def boom(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
            raise RuntimeError(f"{name} (exit 1): crashed")

        monkeypatch.setattr(sy, "run_checked", boom)
        with pytest.raises(RuntimeError, match="sylph sketch failed"):
            sy.identify_by_sylph(["r.fq"], db_dir=str(db))


# ================================================================ CLI


def _k2_stub(top_species: list[dict], calls: list | None = None):
    def identify(files, db_dir=None, threads=4):  # noqa: ANN001 — test stub mirrors real sig
        if calls is not None:
            calls.append((list(files), db_dir, threads))
        return k2.Kraken2Result(
            database={"name": "D5_kraken2", "path": "/db/D5_kraken2"},
            result={
                "top_species": top_species,
                "total_reads": 10000,
                "classified_reads": 100,
                "n_species": len(top_species),
            },
        )

    return identify


def _sy_stub(species: list[dict]):
    def identify(files, db_dir=None, threads=4):  # noqa: ANN001 — test stub mirrors real sig
        return sy.SylphResult(
            database={"name": "D6_sylph", "path": "/db/D6_sylph/db.syldb"},
            result={"species": species, "n_species": len(species)},
        )

    return identify


_K2_SPECIES = [
    {"name": "Escherichia coli", "taxid": 562, "reads": 90, "fraction": 0.5, "confidence": "low"},
    {
        "name": "Salmonella enterica",
        "taxid": 289,
        "reads": 20,
        "fraction": 0.2,
        "confidence": "medium",
    },
]
_SY_SPECIES = [
    {
        "name": "Escherichia coli",
        "ani": 97.5,
        "abundance": 0.80,
        "coverage": 264.0,
        "confidence": "high",
    },
    {
        "name": "Campylobacter jejuni",
        "ani": 92.0,
        "abundance": 0.05,
        "coverage": 10.0,
        "confidence": "low",
    },
]


def _stub_both(monkeypatch, k2_species=None, sy_species=None):
    monkeypatch.setattr(
        k2,
        "identify_by_kraken2",
        _k2_stub(k2_species if k2_species is not None else _K2_SPECIES),
    )
    monkeypatch.setattr(
        sy, "identify_by_sylph", _sy_stub(sy_species if sy_species is not None else _SY_SPECIES)
    )


def _run_capture(monkeypatch, capsys, files, mode="all", **kwargs) -> dict:
    assert _run_meta(list(files), mode, None, **kwargs) == 0
    return json.loads(capsys.readouterr().out)


class TestMetaInputType:
    def test_two_files_detected_as_paired(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["a_1.fq", "a_2.fq"])
        assert payload["input"] == {
            "files": ["a_1.fq", "a_2.fq"],
            "type": "fastq_paired",
        }

    def test_single_fasta_detected_as_contigs(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["contigs.fna"])
        assert payload["input"]["type"] == "fasta"

    def test_single_fastq_detected_as_single(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["reads.fastq"])
        assert payload["input"]["type"] == "fastq_single"

    def test_fa_and_fasta_suffixes_are_contigs(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        for name in ("x.fa", "x.fasta"):
            payload = _run_capture(monkeypatch, capsys, [name])
            assert payload["input"]["type"] == "fasta"


class TestMetaMethods:
    def test_mode_all_runs_both(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["r1.fq", "r2.fq"])
        assert set(payload["methods"]) == {"kraken2", "sylph"}
        assert payload["methods"]["kraken2"]["method"] == "kraken2_bracken"
        assert payload["methods"]["sylph"]["method"] == "sylph_ani"

    def test_mode_kraken2_runs_only_kraken2(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["r.fq"], mode="kraken2")
        assert set(payload["methods"]) == {"kraken2"}
        assert payload["verdict"]["basis"] == ["kraken2"]

    def test_threads_forwarded(self, monkeypatch, capsys):
        calls: list = []
        monkeypatch.setattr(k2, "identify_by_kraken2", _k2_stub(_K2_SPECIES, calls=calls))
        _run_capture(monkeypatch, capsys, ["r.fq"], mode="kraken2", threads=11)
        assert calls == [(["r.fq"], None, 11)]


class TestMetaSylphFastaError:
    def test_contigs_get_structured_sylph_error(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["contigs.fna"])
        assert payload["methods"]["sylph"] == {
            "method": "sylph",
            "error": "sylph requires FASTQ reads, not contigs — use --mode kraken2 for FASTA input",
        }
        assert "result" in payload["methods"]["kraken2"]

    def test_kraken2_only_mode_has_no_sylph_entry(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["contigs.fna"], mode="kraken2")
        assert "sylph" not in payload["methods"]


class TestMetaVerdict:
    def test_merges_and_dedupes_preferring_higher_confidence(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["r1.fq", "r2.fq"])
        top = payload["verdict"]["top_species"]
        names = [row["name"] for row in top]
        assert names == ["Escherichia_coli", "Salmonella_enterica", "Campylobacter_jejuni"]
        assert top[0]["display_name"] == "Escherichia coli"
        e_coli = top[0]
        assert e_coli["method"] == "sylph"  # high (sylph) beats low (kraken2)
        assert e_coli["abundance"] == pytest.approx(0.80)
        assert e_coli["ani"] == 97.5
        assert e_coli["confidence"] == "high"
        assert top[1]["method"] == "kraken2"
        assert top[1]["ani"] is None
        assert payload["verdict"]["richness"] == 3
        assert payload["verdict"]["basis"] == ["kraken2", "sylph"]

    def test_species_group_placeholders_dropped(self, monkeypatch, capsys):
        _stub_both(
            monkeypatch,
            k2_species=[
                {"name": "Campylobacter sp. CCS1377", "fraction": 0.01, "confidence": "low"},
                {"name": "Escherichia coli", "fraction": 0.5, "confidence": "high"},
            ],
            sy_species=[],
        )
        payload = _run_capture(monkeypatch, capsys, ["r.fq"])
        names = [row["name"] for row in payload["verdict"]["top_species"]]
        assert names == ["Escherichia_coli"]
        assert payload["verdict"]["richness"] == 1

    def test_sylph_strain_kept_out_of_the_label(self, monkeypatch, capsys):
        _stub_both(
            monkeypatch,
            k2_species=[],
            sy_species=[
                {
                    "name": "Escherichia coli str. K-12 substr. W3110",
                    "ani": 99.9,
                    "abundance": 0.9,
                    "coverage": 200.0,
                    "confidence": "high",
                }
            ],
        )
        payload = _run_capture(monkeypatch, capsys, ["r.fq"])
        row = payload["verdict"]["top_species"][0]
        assert row["name"] == "Escherichia_coli"
        assert row["display_name"] == "Escherichia coli"
        assert row["strain"] == "str. K-12 substr. W3110"
        assert row["species_raw"] == "Escherichia coli str. K-12 substr. W3110"

    def test_empty_name_entries_skipped(self, monkeypatch, capsys):
        _stub_both(
            monkeypatch,
            k2_species=[{"name": "", "fraction": 0.9, "confidence": "high"}],
        )
        payload = _run_capture(monkeypatch, capsys, ["r.fq"], mode="kraken2")
        assert payload["verdict"]["top_species"] == []
        assert payload["verdict"]["richness"] == 0

    def test_all_methods_failing_still_exits_zero(self, monkeypatch, capsys):
        def boom(files, db_dir=None, threads=4):  # noqa: ANN001 — test stub mirrors real sig
            raise RuntimeError("tool missing")

        monkeypatch.setattr(k2, "identify_by_kraken2", boom)
        monkeypatch.setattr(sy, "identify_by_sylph", boom)
        payload = _run_capture(monkeypatch, capsys, ["r.fq"])
        assert payload["methods"]["kraken2"]["error"] == "tool missing"
        assert payload["methods"]["sylph"]["error"] == "tool missing"
        assert payload["verdict"] == {"top_species": [], "richness": 0, "basis": []}


class TestMetaJsonContract:
    def test_payload_shape(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        payload = _run_capture(monkeypatch, capsys, ["r.fq"])
        assert set(payload) == {
            "analysis_type",
            "tool",
            "version",
            "input",
            "methods",
            "verdict",
        }
        assert payload["analysis_type"] == "metagenomic_profiling"
        assert payload["tool"] == "gside"
        assert payload["version"] == __version__


class TestConfidenceThresholds:
    @pytest.mark.parametrize(
        "fraction,expected",
        [(0.10, "high"), (0.099, "medium"), (0.01, "medium"), (0.009, "low"), (0.0, "low")],
    )
    def test_kraken2_fraction_bands(self, fraction, expected):
        from gside.analysis.kraken2_identifier import _fraction_confidence

        assert _fraction_confidence(fraction) == expected

    @pytest.mark.parametrize(
        "ani,expected",
        [(95.0, "high"), (94.99, "medium"), (93.0, "medium"), (92.99, "low")],
    )
    def test_ani_bands(self, ani, expected):
        from gside.analysis._thresholds import ani_confidence

        assert ani_confidence(ani) == expected


class TestMetaVerdictTieBreak:
    def test_sylph_wins_equal_confidence(self, monkeypatch, capsys):
        _stub_both(
            monkeypatch,
            k2_species=[
                {
                    "name": "Escherichia coli",
                    "taxid": 562,
                    "reads": 100,
                    "fraction": 0.6,
                    "confidence": "high",
                }
            ],
            sy_species=[
                {
                    "name": "Escherichia coli",
                    "ani": 99.0,
                    "abundance": 0.6,
                    "coverage": 100.0,
                    "confidence": "high",
                }
            ],
        )
        payload = _run_capture(monkeypatch, capsys, ["r.fq"])
        row = payload["verdict"]["top_species"][0]
        assert row["method"] == "sylph"
        assert row["ani"] == 99.0


class TestMetaTsv:
    def test_header_and_rows(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        assert _run_meta(["r1.fq", "r2.fq"], "all", None, output_format="tsv") == 0
        lines = capsys.readouterr().out.splitlines()
        assert lines[0] == "method\tname\tabundance\tconfidence\tani"
        assert lines[1] == "sylph\tEscherichia_coli\t0.8\thigh\t97.5"
        assert lines[2] == "kraken2\tSalmonella_enterica\t0.2\tmedium\t"
        assert lines[3] == "sylph\tCampylobacter_jejuni\t0.05\tlow\t92.0"

    def test_kraken2_only_rows_have_blank_ani(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        assert _run_meta(["r.fq"], "kraken2", None, output_format="tsv") == 0
        row = capsys.readouterr().out.splitlines()[1].split("\t")
        assert row[0] == "kraken2" and row[4] == ""

    def test_no_header_omits_header_row(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        assert _run_meta(["r.fq"], "kraken2", None, output_format="tsv", no_header=True) == 0
        lines = capsys.readouterr().out.splitlines()
        assert lines[0].startswith("kraken2\t")
        assert "method\tname" not in lines[0]


class TestMetaOutputFile:
    def test_json_written_to_file_stdout_empty(self, monkeypatch, capsys, tmp_path):
        _stub_both(monkeypatch)
        out = tmp_path / "meta.json"
        assert _run_meta(["r.fq"], "all", None, output=str(out)) == 0
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "Results written to" in captured.err
        assert json.loads(out.read_text())["analysis_type"] == "metagenomic_profiling"

    def test_tsv_written_to_file(self, monkeypatch, capsys, tmp_path):
        _stub_both(monkeypatch)
        out = tmp_path / "meta.tsv"
        assert _run_meta(["r.fq"], "all", None, output_format="tsv", output=str(out)) == 0
        lines = out.read_text().splitlines()
        assert lines[0] == "method\tname\tabundance\tconfidence\tani"
        assert len(lines) == 4


class TestMetaCliCommand:
    def test_registered_and_runs(self, monkeypatch):
        _stub_both(monkeypatch)
        result = CliRunner().invoke(main, ["meta", "r1.fq", "r2.fq"])
        assert result.exit_code == 0
        payload = json.loads(result.output)
        assert payload["input"]["type"] == "fastq_paired"

    def test_help_text(self):
        result = CliRunner().invoke(main, ["meta", "--help"])
        assert result.exit_code == 0
        assert "Metagenomic species profiling" in result.output
        for flag in ("-m", "-D", "-f", "-o", "-t"):
            assert flag in result.output

    def test_invalid_mode_rejected(self):
        result = CliRunner().invoke(main, ["meta", "r.fq", "-m", "bogus"])
        assert result.exit_code == 2

    def test_meta_modes_are_kraken2_sylph_all(self):
        assert cli.META_MODES == ["kraken2", "sylph", "all"]

    def test_default_mode_is_all(self, monkeypatch):
        _stub_both(monkeypatch)
        result = CliRunner().invoke(main, ["meta", "r.fq"])
        payload = json.loads(result.output)
        assert set(payload["methods"]) == {"kraken2", "sylph"}


class TestMetaStructuredErrors:
    def test_missing_db_error_carries_suggested_command(self, monkeypatch, capsys):
        def boom(*args, **kwargs):
            raise RuntimeError("Kraken2 database not found in /x (run: gside db setup --tier meta)")

        monkeypatch.setattr(k2, "identify_by_kraken2", boom)
        payload = _run_capture(monkeypatch, capsys, ["r.fq"], mode="kraken2")
        entry = payload["methods"]["kraken2"]
        assert entry["suggested_command"] == "gside db setup --tier meta"
        assert entry["missing_tier"] == "meta"
        assert entry["help_url"].startswith("https://benlangmead")
        assert entry["estimated_size"] == "~5.5GB"


class TestMetaCsv:
    def test_header_and_rows(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        assert _run_meta(["r1.fq", "r2.fq"], "all", None, output_format="csv") == 0
        lines = capsys.readouterr().out.splitlines()
        assert lines[0] == "method,name,abundance,confidence,ani"
        assert lines[1] == "sylph,Escherichia_coli,0.8,high,97.5"
        assert lines[2] == "kraken2,Salmonella_enterica,0.2,medium,"

    def test_no_header(self, monkeypatch, capsys):
        _stub_both(monkeypatch)
        assert _run_meta(["r.fq"], "kraken2", None, output_format="csv", no_header=True) == 0
        first = capsys.readouterr().out.splitlines()[0]
        assert first.startswith("kraken2,Escherichia_coli,")


class TestSampleManifest:
    def test_header_is_optional(self, tmp_path):
        path = tmp_path / "samples.tsv"
        path.write_text("sample_id\tr1\tr2\nA\tA_1.fq\tA_2.fq\nB\tB.fq\n")
        samples = cli._read_samples(str(path))
        assert [s["id"] for s in samples] == ["A", "B"]
        assert samples[0]["files"] == ["A_1.fq", "A_2.fq"]
        assert samples[1]["files"] == ["B.fq"]

    def test_comments_and_blank_lines_skipped(self, tmp_path):
        path = tmp_path / "samples.tsv"
        path.write_text("# comment\n\nA\ta.fq\n")
        assert [s["id"] for s in cli._read_samples(str(path))] == ["A"]

    def test_duplicate_id_rejected(self, tmp_path):
        path = tmp_path / "samples.tsv"
        path.write_text("A\ta.fq\nA\tb.fq\n")
        with pytest.raises(ValueError, match="duplicate sample id"):
            cli._read_samples(str(path))

    def test_missing_file_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="manifest not found"):
            cli._read_samples(str(tmp_path / "nope.tsv"))

    def test_too_many_columns_rejected(self, tmp_path):
        path = tmp_path / "samples.tsv"
        path.write_text("A\ta.fq\tb.fq\tc.fq\n")
        with pytest.raises(ValueError, match="too many columns"):
            cli._read_samples(str(path))

    def test_batch_preserves_manifest_order_with_jobs(self, monkeypatch, capsys):
        def fake_payload(files, mode, db_dir, threads=4):
            return {"input": {"files": files}, "methods": {}, "verdict": {"top_species": []}}

        monkeypatch.setattr(cli, "_meta_payload", fake_payload)
        samples = [{"id": n, "files": [f"{n}.fq"]} for n in ("a", "b", "c")]
        assert cli._run_meta_samples(samples, "all", None, jobs=3) == 0
        payload = json.loads(capsys.readouterr().out)
        assert [p["sample"] for p in payload] == ["a", "b", "c"]

    def test_batch_tsv_prefixes_sample_column(self, monkeypatch, capsys):
        def fake_payload(files, mode, db_dir, threads=4):
            return {
                "methods": {},
                "verdict": {
                    "top_species": [
                        {
                            "method": "sylph",
                            "name": "Escherichia_coli",
                            "abundance": 0.5,
                            "confidence": "high",
                            "ani": 99.0,
                        }
                    ]
                },
            }

        monkeypatch.setattr(cli, "_meta_payload", fake_payload)
        samples = [{"id": "A", "files": ["a.fq"]}, {"id": "B", "files": ["b.fq"]}]
        assert cli._run_meta_samples(samples, "sylph", None, output_format="tsv") == 0
        lines = capsys.readouterr().out.splitlines()
        assert lines[0] == "sample\tmethod\tname\tabundance\tconfidence\tani"
        assert lines[1] == "A\tsylph\tEscherichia_coli\t0.5\thigh\t99.0"
        assert lines[2].startswith("B\tsylph\t")

    def test_jobs_without_samples_is_a_usage_error(self):
        result = CliRunner().invoke(main, ["meta", "r.fq", "-j", "2"])
        assert result.exit_code == 2
        assert "--jobs only applies" in result.output

    def test_reads_and_samples_are_exclusive(self, tmp_path):
        path = tmp_path / "s.tsv"
        path.write_text("A\ta.fq\n")
        result = CliRunner().invoke(main, ["meta", "r.fq", "--samples", str(path)])
        assert result.exit_code == 2
        assert "not both" in result.output

    def test_missing_input_is_a_usage_error(self):
        result = CliRunner().invoke(main, ["meta"])
        assert result.exit_code == 2
        assert "missing input" in result.output


class TestSampleManifestCliDispatch:
    def test_cli_runs_manifest_samples(self, monkeypatch, tmp_path, capsys):
        def fake_payload(files, mode, db_dir, threads=4):
            return {"input": {"files": files}, "methods": {}, "verdict": {"top_species": []}}

        monkeypatch.setattr(cli, "_meta_payload", fake_payload)
        manifest = tmp_path / "samples.tsv"
        manifest.write_text("sample_id\tr1\tr2\nA\ta_1.fq\ta_2.fq\nB\tb.fq\n")
        result = CliRunner().invoke(main, ["meta", "--samples", str(manifest), "-j", "2"])
        assert result.exit_code == 0, result.output
        payload = json.loads(result.output)
        assert [p["sample"] for p in payload] == ["A", "B"]

    def test_row_missing_reads_column_rejected(self, tmp_path):
        path = tmp_path / "samples.tsv"
        path.write_text("A\n")
        with pytest.raises(ValueError, match="expected sample_id"):
            cli._read_samples(str(path))

    def test_empty_manifest_rejected(self, tmp_path):
        path = tmp_path / "samples.tsv"
        path.write_text("# only comments\n\n")
        with pytest.raises(ValueError, match="no samples found"):
            cli._read_samples(str(path))

    def test_cli_rejects_invalid_manifest(self, tmp_path):
        path = tmp_path / "samples.tsv"
        path.write_text("A\ta.fq\nA\tb.fq\n")
        result = CliRunner().invoke(main, ["meta", "--samples", str(path)])
        assert result.exit_code == 2
        assert "duplicate sample id" in result.output
