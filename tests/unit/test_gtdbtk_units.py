"""Pure unit tests for gtdbtk_identifier (subprocess stubbed, no binaries)."""

from __future__ import annotations

import types
from pathlib import Path

import pytest

import gside.analysis.gtdbtk_identifier as gt
from gside.analysis.gtdbtk_identifier import (
    _species_from_classification,
    identify_by_gtdbtk,
)

_BAC120_LINE = (
    "d__Bacteria;p__Firmicutes;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;"
    "g__Clostridium;s__Clostridium perfringens"
)


def _fake_run_checked(
    summary: str = "gtdbtk.bac120.summary.tsv",
    classification: str = _BAC120_LINE,
    cmds: list | None = None,
    genomes: dict[str, str] | None = None,
):
    """Fake run_checked mimicking gtdbtk: writes the summary under --out_dir."""

    def run_checked(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
        if cmds is not None:
            cmds.append((list(cmd), timeout, name, kwargs))
        genome_dir = Path(cmd[cmd.index("--genome_dir") + 1])
        if genomes is not None:
            genomes.update({p.name: p.read_text() for p in genome_dir.iterdir()})
        out_dir = Path(cmd[cmd.index("--out_dir") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / summary).write_text(f"user_genome\tclassification\nquery\t{classification}\n")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    return run_checked


def _ready(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(gt, "GTDB_DB", tmp_path / "gtdb")
    monkeypatch.setattr(gt, "which", lambda name: "/bin/gtdbtk")


def _contigs(tmp_path: Path) -> Path:
    q = tmp_path / "q.fna"
    q.write_text(">q\nACGTACGT\n")
    return q


class TestToDict:
    def test_shape(self):
        d = gt.GtdbkResult(database={"name": "GTDB-Tk"}, result={"species": "X"}).to_dict()
        assert d["analysis_type"] == "species_identification"
        assert d["method"] == "gtdbtk"
        assert d["database"] == {"name": "GTDB-Tk"}
        assert d["result"] == {"species": "X"}


class TestSpeciesParsing:
    def test_extracts_last_s_rank(self):
        line = "d__Bacteria;p__Firmicutes;s__Clostridium perfringens"
        assert _species_from_classification(line) == "Clostridium perfringens"

    def test_no_s_rank_is_unknown(self):
        assert _species_from_classification("d__Bacteria;p__Firmicutes;g__Clostridium") == "Unknown"

    def test_empty_s_rank_is_unknown(self):
        assert _species_from_classification("d__Bacteria;s__") == "Unknown"

    def test_empty_string_is_unknown(self):
        assert _species_from_classification("") == "Unknown"

    def test_tolerates_rank_whitespace(self):
        assert _species_from_classification("d__Bacteria; s__Vibrio parahaemolyticus") == (
            "Vibrio parahaemolyticus"
        )


class TestUnavailable:
    def test_missing_db_raises_actionable(self, monkeypatch):
        monkeypatch.setattr(gt, "GTDB_DB", None)
        monkeypatch.setattr(gt, "which", lambda name: "/bin/gtdbtk")
        with pytest.raises(RuntimeError, match="GTDB-Tk not available") as exc:
            identify_by_gtdbtk("q.fna")
        assert "GTDBTK_DATA_PATH" in str(exc.value)
        assert "https://ecogenomics.github.io/GTDBTk/" in str(exc.value)

    def test_missing_binary_raises_actionable(self, monkeypatch, tmp_path):
        monkeypatch.setattr(gt, "GTDB_DB", tmp_path)
        monkeypatch.setattr(gt, "which", lambda name: None)
        with pytest.raises(RuntimeError, match="GTDB-Tk not available") as exc:
            identify_by_gtdbtk("q.fna")
        assert "GTDBTK_DATA_PATH" in str(exc.value)


class TestClassifyWf:
    def test_bac120_happy_path(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)
        monkeypatch.setattr(gt, "run_checked", _fake_run_checked())
        r = identify_by_gtdbtk(_contigs(tmp_path))
        assert r.method == "gtdbtk"
        assert r.result["species"] == "Clostridium perfringens"
        assert r.result["confidence"] == "high"
        assert r.result["classification"] == _BAC120_LINE

    def test_ar122_fallback(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)
        monkeypatch.setattr(
            gt,
            "run_checked",
            _fake_run_checked(summary="gtdbtk.ar122.summary.tsv", classification="d__Archaea;s__X"),
        )
        r = identify_by_gtdbtk(_contigs(tmp_path))
        assert r.result["species"] == "X"
        assert r.result["classification"] == "d__Archaea;s__X"

    def test_command_shape_genome_copy_and_db(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)
        cmds: list = []
        genomes: dict[str, str] = {}
        monkeypatch.setattr(gt, "run_checked", _fake_run_checked(cmds=cmds, genomes=genomes))
        identify_by_gtdbtk(_contigs(tmp_path), threads=9)
        cmd, timeout, name, _kwargs = cmds[0]
        assert cmd[0] == "/bin/gtdbtk" and cmd[1] == "classify_wf"
        assert cmd[cmd.index("--database") + 1] == str(tmp_path / "gtdb")
        assert cmd[cmd.index("--cpus") + 1] == "9"
        assert timeout == 3600
        assert name == "gtdbtk"
        assert genomes == {"query.fna": ">q\nACGTACGT\n"}

    def test_genome_dir_is_separate_from_out_dir(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)
        seen: list[tuple[Path, Path]] = []

        def run_checked(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
            genome_dir = Path(cmd[cmd.index("--genome_dir") + 1])
            out_dir = Path(cmd[cmd.index("--out_dir") + 1])
            seen.append((genome_dir, out_dir))
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "gtdbtk.bac120.summary.tsv").write_text(
                f"user_genome\tclassification\nq\t{_BAC120_LINE}\n"
            )
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(gt, "run_checked", run_checked)
        identify_by_gtdbtk(_contigs(tmp_path))
        genome_dir, out_dir = seen[0]
        assert genome_dir != out_dir and out_dir.name == "output"

    def test_failure_propagates(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)

        def boom(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
            raise RuntimeError("gtdbtk (exit 1): crashed")

        monkeypatch.setattr(gt, "run_checked", boom)
        with pytest.raises(RuntimeError, match="gtdbtk"):
            identify_by_gtdbtk(_contigs(tmp_path))

    def test_no_summary_is_unknown(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)

        def run_checked(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(gt, "run_checked", run_checked)
        r = identify_by_gtdbtk(_contigs(tmp_path))
        assert r.result["species"] == "Unknown"
        assert r.result["confidence"] == "low"
        assert r.result["classification"] == ""

    def test_header_only_summary_is_unknown(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)

        def run_checked(cmd: list[str], *, timeout: int, name: str, **kwargs: object):
            out_dir = Path(cmd[cmd.index("--out_dir") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "gtdbtk.bac120.summary.tsv").write_text("user_genome\tclassification\n")
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(gt, "run_checked", run_checked)
        r = identify_by_gtdbtk(_contigs(tmp_path))
        assert r.result["species"] == "Unknown"
        assert r.result["classification"] == ""

    def test_classification_without_species_rank(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)
        monkeypatch.setattr(
            gt,
            "run_checked",
            _fake_run_checked(classification="d__Bacteria;p__Firmicutes;g__Clostridium"),
        )
        r = identify_by_gtdbtk(_contigs(tmp_path))
        assert r.result["species"] == "Unknown"
        assert r.result["confidence"] == "low"
        assert r.result["classification"] == "d__Bacteria;p__Firmicutes;g__Clostridium"


class TestTempDirIntermediates:
    def test_genome_copied_not_moved(self, monkeypatch, tmp_path):
        _ready(monkeypatch, tmp_path)
        q = _contigs(tmp_path)
        monkeypatch.setattr(gt, "run_checked", _fake_run_checked())
        identify_by_gtdbtk(q)
        assert q.read_text() == ">q\nACGTACGT\n"
