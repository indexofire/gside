"""Offline unit tests for gside/db.py (tmp_path + monkeypatch, no network)."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import shutil
import tarfile
from pathlib import Path

import pytest

from gside import db
from gside.db import (
    _build_panel_from_manifest,
    _check_mash,
    _check_panel,
    _copy_from_source,
    _download_mash_zenodo,
    _try_download_panel_release,
    db_setup,
    db_status,
    run_db_command,
)


@pytest.fixture()
def root(monkeypatch, tmp_path):
    monkeypatch.setattr(db, "SPECIES_DB_DIR", tmp_path / "db")
    return tmp_path / "db"


def _touch(path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(b"x")
    return p


class TestStatus:
    def test_empty_all_not_ready(self, root):
        status = db_status()
        assert set(status) == {
            "D1_marker",
            "D2_ani",
            "D3_mash",
            "D4_sourmash",
        }
        assert [status[k]["tier"] for k in ("D1_marker", "D2_ani", "D3_mash")] == [
            "mini",
            "panel",
            "mash",
        ]
        assert not any(v["ready"] for v in status.values())

    def test_ready_branches(self, root):
        _touch(root / "D1_marker" / "markers.fasta")
        _touch(root / "D2_ani" / "panel.sketch" / "sketches.db")
        _touch(root / "D3_mash" / "payload.bin")
        _touch(root / "D4_sourmash" / "gtdb-reps-k31.zip")
        _touch(root / "D4_sourmash" / "lineages.csv")
        status = db_status()
        assert all(v["ready"] for v in status.values())

    def test_panel_sketch_as_file(self, root):
        _touch(root / "D2_ani" / "panel.sketch")
        assert _check_panel()

    def test_mash_msh(self, root):
        _touch(root / "D3_mash" / "mash.msh")
        assert _check_mash()


class TestSetupTiers:
    def test_mini_ships_with_repo(self):
        assert db_setup("mini") == {"info": "mini tier ships with repo"}

    def test_unknown_tier(self):
        with pytest.raises(ValueError, match="unknown tier"):
            db_setup("wat")

    def test_already_ready(self, root):
        _touch(root / "D2_ani" / "panel.sketch" / "sketches.db")
        _touch(root / "D3_mash" / "mash.msh")
        assert db_setup("panel") == {"D2_ani": "already ready"}
        assert db_setup("mash") == {"D3_mash": "already ready"}

    def test_all_tier_covers_panel_and_mash(self, root):
        _touch(root / "D2_ani" / "panel.sketch" / "sketches.db")
        _touch(root / "D3_mash" / "mash.msh")
        assert db_setup("all") == {"D2_ani": "already ready", "D3_mash": "already ready"}


class TestCopyFromSource:
    def test_copy_exact_dir(self, root, tmp_path):
        src = tmp_path / "D3_mash"
        _touch(src / "mash.msh")
        assert _copy_from_source(str(src), "D3_mash") == f"copied from {src}"
        assert (root / "D3_mash" / "mash.msh").is_file()

    def test_copy_parent_appends_name(self, root, tmp_path):
        _touch(tmp_path / "D2_ani" / "panel.sketch" / "sketches.db")
        assert _copy_from_source(str(tmp_path), "D2_ani").startswith("copied from")
        assert (root / "D2_ani" / "panel.sketch" / "sketches.db").is_file()

    def test_missing_source(self, root, tmp_path):
        msg = _copy_from_source(str(tmp_path / "nope"), "D2_ani")
        assert msg.startswith("ERROR: source not found")


def _write_gz(path: Path, payload: bytes = b"fake-sketch") -> None:
    with gzip.open(path, "wb") as fh:
        fh.write(payload)


class TestMashDownload:
    def test_md5_mismatch_cleans_up(self, root, monkeypatch, tmp_path):
        dst = root / "D3_mash"
        dst.mkdir(parents=True)

        def fake_download(url, dst):
            _write_gz(Path(dst), b"junk")

        monkeypatch.setattr(db, "_download_file", fake_download)
        msg = _download_mash_zenodo(dst)
        assert msg.startswith("ERROR: MD5 mismatch")
        assert not (root / "D3_mash" / "mash.msh.gz").exists()

    def test_happy_path(self, root, monkeypatch, tmp_path):
        dst = root / "D3_mash"
        dst.mkdir(parents=True)
        buf = io.BytesIO()
        with gzip.open(buf, "wb") as fh:
            fh.write(b"fake-sketch-bytes")
        gz_bytes = buf.getvalue()

        def fake_download(url, dst):
            Path(dst).write_bytes(gz_bytes)

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "MASH_MD5", hashlib.md5(gz_bytes).hexdigest())
        msg = _download_mash_zenodo(dst)
        assert msg == "downloaded from Zenodo (community sketch)"
        assert (dst / "mash.msh").read_bytes() == b"fake-sketch-bytes"

    def test_network_failure(self, root, monkeypatch):
        def boom(url, dst):
            raise OSError("offline")

        monkeypatch.setattr(db, "_download_file", boom)
        assert _download_mash_zenodo(root / "D3_mash").startswith("ERROR: Zenodo")


def _write_tar(path: Path, names: list[str]) -> None:
    with tarfile.open(path, "w:gz") as tf:
        for name in names:
            info = tarfile.TarInfo(name)
            info.size = 1
            tf.addfile(info, io.BytesIO(b"x"))


class TestPanelRelease:
    def test_valid_tarball(self, root, monkeypatch, tmp_path):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)

        def fake_download(url, dst):
            _write_tar(Path(dst), ["panel.sketch/sketches.db"])

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "PANEL_SHA256", "")
        msg = _try_download_panel_release(dst)
        assert msg == "downloaded from GitHub Release (pre-built sketch)"

    def test_invalid_tarball(self, root, monkeypatch):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)

        def fake_download(url, dst):
            _write_tar(Path(dst), ["other/file.txt"])

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "PANEL_SHA256", "")
        msg = _try_download_panel_release(dst)
        assert msg.startswith("ERROR: downloaded archive did not contain")

    def test_download_failure_returns_none(self, root, monkeypatch):
        def boom(url, dst):
            raise OSError("offline")

        monkeypatch.setattr(db, "_download_file", boom)
        assert _try_download_panel_release(root / "D2_ani") is None


class TestBuildFromManifest:
    def test_missing_manifest(self, monkeypatch, tmp_path):
        monkeypatch.setattr(db, "_MANIFEST_DIR", tmp_path / "empty")
        assert _build_panel_from_manifest(tmp_path / "dst").startswith(
            "ERROR: panel manifest not found"
        )

    def test_missing_datasets_cli(self, monkeypatch, tmp_path):
        manifest = tmp_path / "manifests"
        manifest.mkdir()
        (manifest / "panel_accessions.tsv").write_text("accession\n")
        (manifest / "metadata.tsv").write_text("x\n")
        monkeypatch.setattr(db, "_MANIFEST_DIR", manifest)
        monkeypatch.setattr(db, "which", lambda tool: None)
        assert _build_panel_from_manifest(tmp_path / "dst").startswith(
            "ERROR: NCBI datasets CLI not found"
        )

    def test_missing_skani(self, monkeypatch, tmp_path):
        manifest = tmp_path / "manifests"
        manifest.mkdir()
        (manifest / "panel_accessions.tsv").write_text("accession\n")
        (manifest / "metadata.tsv").write_text("x\n")
        monkeypatch.setattr(db, "_MANIFEST_DIR", manifest)
        monkeypatch.setattr(
            db, "which", lambda tool: "/bin/datasets" if tool == "datasets" else None
        )
        assert _build_panel_from_manifest(tmp_path / "dst").startswith("ERROR: skani not found")


class TestFormatStatus:
    def test_table_matches_print_status_bytes(self, root, capsys):
        _touch(root / "D1_marker" / "markers.fasta")
        run_db_command(["status"])
        expected = capsys.readouterr().out
        assert db.format_status("table") + "\n" == expected

    def test_tsv_rows(self, root):
        _touch(root / "D1_marker" / "markers.fasta")
        text = db.format_status("tsv")
        lines = text.splitlines()
        assert lines[0] == "COMPONENT\tREADY\tTIER\tINPUTS\tPATH"
        assert lines[1].startswith("D1_marker\tready\tmini\tFASTA/FASTQ\t")
        assert lines[2].split("\t")[:4] == ["D2_ani", "missing", "panel", "FASTA"]
        assert "✅" not in text
        assert "Run 'gside" not in text

    def test_json_roundtrip(self, root):
        assert json.loads(db.format_status("json")) == db_status()

    def test_unknown_format_raises(self):
        with pytest.raises(ValueError, match="unknown format"):
            db.format_status("wat")


class TestFormatTiers:
    def test_table_matches_run_db_command_bytes(self, capsys):
        run_db_command(["list"])
        expected = capsys.readouterr().out
        assert db.format_tiers("table") + "\n" == expected

    def test_table_has_bordered_header_and_all_tiers(self):
        text = db.format_tiers("table")
        assert "gside database tiers" in text
        assert "│ TIER" in text
        assert "│ READY" in text
        assert "│ SIZE" in text
        assert "│ INPUTS" in text
        assert "│ DESCRIPTION" in text
        assert "┌" in text and "└" in text
        for name in db._TIERS:
            assert f"│ {name}" in text
        ready_values = {
            line.split("│")[2].strip()
            for line in text.splitlines()
            if line.startswith("│ ") and line.split("│")[1].strip()
        }
        assert ready_values <= {"ready", "missing", "READY"}

    def test_tsv_rows(self):
        lines = db.format_tiers("tsv").splitlines()
        assert lines[0] == "NAME\tREADY\tSIZE\tINPUTS\tDESCRIPTION"
        assert len(lines) == 1 + len(db._TIERS)
        assert lines[1].split("\t")[:4] == ["mini", "ready", "~1MB", "FASTA/FASTQ"]

    def test_tsv_ready_aggregates_tier_items(self, root):
        _touch(root / "D2_ani" / "panel.sketch" / "sketches.db")
        lines = db.format_tiers("tsv").splitlines()
        by_name = {ln.split("\t")[0]: ln.split("\t")[1] for ln in lines[1:]}
        assert by_name["panel"] == "ready"
        assert by_name["all"] == "missing"
        assert by_name["mash"] == "missing"

    def test_json_roundtrip_with_ready(self):
        data = json.loads(db.format_tiers("json"))
        assert set(data) == {*db._TIERS, "custom_databases"}
        assert data["mini"]["ready"] is True
        assert isinstance(data["panel"]["ready"], bool)

    def test_csv_quotes_commas(self):
        lines = db.format_tiers("csv").splitlines()
        assert lines[0] == "NAME,READY,SIZE,INPUTS,DESCRIPTION"
        sourmash_row = next(ln for ln in lines[1:] if ln.startswith("sourmash,"))
        assert '"GTDB gather database (sourmash, sketch k=31)"' in sourmash_row

    def test_unknown_format_raises(self):
        with pytest.raises(ValueError, match="unknown format"):
            db.format_tiers("wat")


class TestRunDbCommand:
    def test_status_and_list(self, root, capsys):
        assert run_db_command(["status"]) == 0
        assert run_db_command(["list"]) == 0
        assert "D3_mash" in capsys.readouterr().out

    def test_unknown_subcommand(self):
        assert run_db_command(["wat"]) == 1

    def test_setup_error_propagates_exit_1(self, monkeypatch):
        monkeypatch.setattr(db, "db_setup", lambda tier, source: {"x": "ERROR: boom"})
        assert run_db_command(["setup", "panel"]) == 1

    def test_setup_ok(self, root):
        _touch(root / "D2_ani" / "panel.sketch" / "sketches.db")
        assert run_db_command(["setup", "panel"]) == 0


class TestSourmashTier:
    def test_already_ready(self, root):
        _touch(root / "D4_sourmash" / "gtdb-reps-k31.zip")
        _touch(root / "D4_sourmash" / "lineages.csv")
        assert db_setup("sourmash") == {"D4_sourmash": "already ready"}

    def test_status_row(self, root):
        _touch(root / "D4_sourmash" / "gtdb-reps-k31.zip")
        _touch(root / "D4_sourmash" / "lineages.csv")
        assert db_status()["D4_sourmash"]["ready"]

    def test_download_failure(self, root, monkeypatch):
        def boom(url, dst):
            raise OSError("offline")

        monkeypatch.setattr(db, "_download_file", boom)
        msg = db_setup("sourmash")["D4_sourmash"]
        assert msg.startswith("ERROR: farm download failed")

    def test_download_renames_in_place(self, root, monkeypatch, tmp_path):
        def fake_download(url, dst):
            Path(dst).write_bytes(b"x")

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "SOURMASH_SIG_SHA256", hashlib.sha256(b"x").hexdigest())
        monkeypatch.setattr(db, "SOURMASH_LINEAGES_SHA256", hashlib.sha256(b"x").hexdigest())
        msg = db_setup("sourmash")["D4_sourmash"]
        assert msg.startswith("downloaded from farm.cse.ucdavis.edu")
        assert (root / "D4_sourmash" / "gtdb-reps-k31.zip").is_file()
        assert (root / "D4_sourmash" / "lineages.csv").is_file()

    def test_copy_from_source(self, root, tmp_path):
        src = tmp_path / "D4_sourmash"
        _touch(src / "gtdb-reps-k31.zip")
        _touch(src / "lineages.csv")
        assert db_setup("sourmash", source=str(tmp_path))["D4_sourmash"].startswith("copied from")


class TestSetupDownloadPaths:
    def test_panel_download_then_build(self, root, monkeypatch):
        calls = []
        monkeypatch.setattr(db, "_try_download_panel_release", lambda dst: None)
        monkeypatch.setattr(
            db, "_build_panel_from_manifest", lambda dst: calls.append(dst) or "built!"
        )
        assert db_setup("panel")["D2_ani"] == "built!"
        assert calls[0] == root / "D2_ani"

    def test_panel_download_success(self, root, monkeypatch):
        monkeypatch.setattr(db, "_try_download_panel_release", lambda dst: "downloaded!")
        assert db_setup("panel")["D2_ani"] == "downloaded!"

    def test_panel_release_publishes_taxa_map(self, root, monkeypatch, tmp_path):
        manifests = tmp_path / "manifests"
        manifests.mkdir()
        (manifests / "metadata.tsv").write_text("fna\tspecies\nGCF_1.fna\tSpecies one\n")
        monkeypatch.setattr(db, "_MANIFEST_DIR", manifests)
        monkeypatch.setattr(db, "_try_download_panel_release", lambda dst: "downloaded!")
        db_setup("panel")
        assert (root / "D2_ani" / "metadata.tsv").is_file()

    def test_mash_download_path(self, root, monkeypatch):
        monkeypatch.setattr(db, "_download_mash_zenodo", lambda dst: "dl!")
        assert db_setup("mash")["D3_mash"] == "dl!"

    def test_sourmash_download_path(self, root, monkeypatch):
        monkeypatch.setattr(db, "_download_sourmash_farm", lambda dst: "dl!")
        assert db_setup("sourmash")["D4_sourmash"] == "dl!"


class TestRunSetupArgs:
    def test_tier_and_source(self, monkeypatch):
        seen = {}
        monkeypatch.setattr(
            db, "db_setup", lambda tier, source="": seen.update(tier=tier, source=source) or {}
        )
        assert run_db_command(["setup", "--tier", "mash", "--source", "/s"]) == 0
        assert seen == {"tier": "mash", "source": "/s"}

    def test_bare_setup_defaults_panel(self, monkeypatch):
        seen = {}
        monkeypatch.setattr(db, "db_setup", lambda tier, source="": seen.update(tier=tier) or {})
        assert run_db_command(["setup"]) == 0
        assert seen == {"tier": "panel"}

    def test_unknown_args_skipped(self, monkeypatch):
        monkeypatch.setattr(db, "db_setup", lambda tier, source="": {})
        assert run_db_command(["setup", "--bogus", "x", "--tier"]) == 0


class TestSetupSourcePaths:
    def test_panel_source(self, root, tmp_path):
        _touch(tmp_path / "src" / "D2_ani" / "panel.sketch" / "sketches.db")
        assert db_setup("panel", source=str(tmp_path / "src"))["D2_ani"].startswith("copied from")

    def test_mash_source(self, root, tmp_path):
        _touch(tmp_path / "src" / "D3_mash" / "mash.msh")
        assert db_setup("mash", source=str(tmp_path / "src"))["D3_mash"].startswith("copied from")


class TestSourmashRenameErrors:
    def test_rename_failure(self, root, monkeypatch):
        def fake_download(url, dst):
            Path(dst).write_bytes(b"x")

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "SOURMASH_SIG_SHA256", hashlib.sha256(b"x").hexdigest())
        monkeypatch.setattr(db, "SOURMASH_LINEAGES_SHA256", hashlib.sha256(b"x").hexdigest())
        dst = root / "D4_sourmash"
        dst.mkdir(parents=True)
        # Pre-existing directory at the rename target makes os.rename fail.
        (dst / "lineages.csv").mkdir()
        msg = db_setup("sourmash")["D4_sourmash"]
        assert msg.startswith("ERROR:")

    def test_incomplete(self, root, monkeypatch):
        def fake_download(url, dst):
            Path(dst).write_bytes(b"x")

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "SOURMASH_SIG_SHA256", hashlib.sha256(b"x").hexdigest())
        monkeypatch.setattr(db, "SOURMASH_LINEAGES_SHA256", hashlib.sha256(b"x").hexdigest())
        monkeypatch.setattr(db, "_check_sourmash", lambda: False)
        assert db_setup("sourmash")["D4_sourmash"].startswith("ERROR: download incomplete")


class TestPinnedChecksumConstants:
    def test_constants(self):
        assert db.PANEL_SHA256 == (
            "f9bdbeeaa744a6ed9fafef98e854d309eebc8d060df787345c2d7cdd1e36e595"
        )
        assert db.SOURMASH_SIG_SHA256 == (
            "6f3cc74f0bdcf6f84466ba91f87d6f2bf30081d5814bdf006cafd1e9cbbeada9"
        )
        assert db.SOURMASH_LINEAGES_SHA256 == (
            "98bceab27a50f08b2f777ca7bdfb57c88aabe5ce1fa54cfb103dd0ad51b67624"
        )


class TestVerifySha256:
    def test_matching_digest_passes(self, tmp_path):
        artifact = tmp_path / "artifact.bin"
        artifact.write_bytes(b"payload")
        expected = hashlib.sha256(b"payload").hexdigest()
        assert db._verify_sha256(artifact, expected) is True

    def test_multichunk_streaming_matches(self, tmp_path):
        payload = bytes(range(256)) * 12288
        artifact = tmp_path / "big.bin"
        artifact.write_bytes(payload)
        expected = hashlib.sha256(payload).hexdigest()
        assert db._verify_sha256(artifact, expected) is True

    def test_wrong_digest_fails(self, tmp_path):
        artifact = tmp_path / "artifact.bin"
        artifact.write_bytes(b"payload")
        assert db._verify_sha256(artifact, "0" * 64) is False

    def test_empty_expected_warns_and_passes(self, tmp_path, capsys):
        artifact = tmp_path / "artifact.bin"
        artifact.write_bytes(b"payload")
        assert db._verify_sha256(artifact, "") is True
        out = capsys.readouterr().out
        assert "WARNING: no pinned checksum for artifact.bin" in out
        assert "skipping integrity verification" in out

    def test_missing_file_fails_closed(self, tmp_path):
        assert db._verify_sha256(tmp_path / "missing.bin", "0" * 64) is False


class TestTarExtractionSafety:
    def test_traversal_member_rejected(self, root, monkeypatch):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)

        def fake_download(url, d):
            _write_tar(Path(d), ["panel.sketch/sketches.db", "../evil.txt"])

        monkeypatch.setattr(db, "_download_file", fake_download)
        msg = _try_download_panel_release(dst)
        assert isinstance(msg, str) and msg.startswith("ERROR:")
        assert not (root / "evil.txt").exists()

    def test_absolute_path_member_rejected(self, root, monkeypatch):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)

        def fake_download(url, d):
            _write_tar(Path(d), ["panel.sketch/sketches.db", "/../evil-abs.txt"])

        monkeypatch.setattr(db, "_download_file", fake_download)
        msg = _try_download_panel_release(dst)
        assert isinstance(msg, str) and msg.startswith("ERROR:")
        assert not (root / "evil-abs.txt").exists()

    def test_benign_nested_tar_still_extracts(self, root, monkeypatch):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)

        def fake_download(url, d):
            _write_tar(Path(d), ["panel.sketch/sketches.db", "panel.sketch/nested/deep/file.txt"])

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "PANEL_SHA256", "")
        msg = _try_download_panel_release(dst)
        assert msg == "downloaded from GitHub Release (pre-built sketch)"
        assert (dst / "panel.sketch" / "sketches.db").is_file()
        assert (dst / "panel.sketch" / "nested" / "deep" / "file.txt").is_file()


def _tar_bytes(names: list[str]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for name in names:
            info = tarfile.TarInfo(name)
            info.size = 1
            tf.addfile(info, io.BytesIO(b"x"))
    return buf.getvalue()


class TestPanelChecksumVerification:
    def test_pinned_match_extracts(self, root, monkeypatch):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)
        archive_bytes = _tar_bytes(["panel.sketch/sketches.db"])

        def fake_download(url, d):
            Path(d).write_bytes(archive_bytes)

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "PANEL_SHA256", hashlib.sha256(archive_bytes).hexdigest())
        msg = _try_download_panel_release(dst)
        assert msg == "downloaded from GitHub Release (pre-built sketch)"
        assert (dst / "panel.sketch" / "sketches.db").is_file()

    def test_pinned_mismatch_fails_closed(self, root, monkeypatch):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)

        def fake_download(url, d):
            _write_tar(Path(d), ["panel.sketch/sketches.db"])

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "PANEL_SHA256", "f" * 64)
        msg = _try_download_panel_release(dst)
        assert msg is not None and msg.startswith("ERROR: SHA256 mismatch")
        assert not (dst / "panel.sketch.tar.gz").exists()
        assert not (dst / "panel.sketch").exists()

    def test_unpinned_warns_and_continues(self, root, monkeypatch, capsys):
        dst = root / "D2_ani"
        dst.mkdir(parents=True)

        def fake_download(url, d):
            _write_tar(Path(d), ["panel.sketch/sketches.db"])

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "PANEL_SHA256", "")
        msg = _try_download_panel_release(dst)
        assert msg == "downloaded from GitHub Release (pre-built sketch)"
        assert "WARNING: no pinned checksum for panel.sketch.tar.gz" in capsys.readouterr().out


class TestSourmashChecksumVerification:
    def test_lineages_mismatch_fails_closed(self, root, monkeypatch):
        def fake_download(url, d):
            Path(d).write_bytes(b"x")

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "SOURMASH_SIG_SHA256", hashlib.sha256(b"x").hexdigest())
        dst = root / "D4_sourmash"
        dst.mkdir(parents=True)
        msg = db._download_sourmash_farm(dst)
        assert msg.startswith("ERROR: SHA256 mismatch for gtdb-rs226-reps.lineages.csv")
        assert not (dst / "gtdb-rs226-reps.lineages.csv").exists()
        assert not (dst / "lineages.csv").exists()
        assert not (dst / "gtdb-reps-k31.zip").exists()

    def test_sig_unpinned_warns_and_lineages_verified(self, root, monkeypatch, capsys):
        def fake_download(url, d):
            Path(d).write_bytes(b"x")

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "SOURMASH_SIG_SHA256", "")
        monkeypatch.delenv("GSIDE_SOURMASH_SIG_SHA256", raising=False)
        monkeypatch.setattr(db, "SOURMASH_LINEAGES_SHA256", hashlib.sha256(b"x").hexdigest())
        dst = root / "D4_sourmash"
        dst.mkdir(parents=True)
        msg = db._download_sourmash_farm(dst)
        assert msg.startswith("downloaded from farm.cse.ucdavis.edu")
        out = capsys.readouterr().out
        assert "WARNING: no pinned checksum for gtdb-rs226-reps.k31.sig.zip" in out
        assert "WARNING: no pinned checksum for gtdb-rs226-reps.lineages.csv" not in out


class TestMashStreamingMd5:
    def test_multichunk_hash_without_whole_file_read(self, root, monkeypatch):
        dst = root / "D3_mash"
        dst.mkdir(parents=True)
        payload = b"a" * (2 * 1024 * 1024 + 123)
        buf = io.BytesIO()
        with gzip.open(buf, "wb") as fh:
            fh.write(payload)
        gz_bytes = buf.getvalue()

        def fake_download(url, d):
            Path(d).write_bytes(gz_bytes)

        def no_read_bytes(self):
            raise AssertionError("whole-file read_bytes; MD5 must stream in chunks")

        monkeypatch.setattr(db, "_download_file", fake_download)
        monkeypatch.setattr(db, "MASH_MD5", hashlib.md5(gz_bytes).hexdigest())
        monkeypatch.setattr(Path, "read_bytes", no_read_bytes)
        msg = _download_mash_zenodo(dst)
        monkeypatch.undo()
        assert msg == "downloaded from Zenodo (community sketch)"
        assert (dst / "mash.msh").read_bytes() == payload
        assert not (dst / "mash.msh.gz").exists()


class TestSourmashChecksumOverride:
    def test_env_override_rejects_wrong_digest(self, root, monkeypatch):
        monkeypatch.setattr(
            db, "_download_file", lambda url, dst: Path(dst).write_bytes(b"payload")
        )
        monkeypatch.setenv("GSIDE_SOURMASH_SIG_SHA256", "0" * 64)
        msg = db_setup("sourmash")["D4_sourmash"]
        assert msg.startswith("ERROR: SHA256 mismatch")

    def test_env_override_accepts_matching_digest(self, root, monkeypatch):
        monkeypatch.setattr(
            db, "_download_file", lambda url, dst: Path(dst).write_bytes(b"payload")
        )
        import hashlib as _hashlib

        monkeypatch.setenv("GSIDE_SOURMASH_SIG_SHA256", _hashlib.sha256(b"payload").hexdigest())
        monkeypatch.setattr(db, "SOURMASH_LINEAGES_SHA256", _hashlib.sha256(b"payload").hexdigest())
        msg = db_setup("sourmash")["D4_sourmash"]
        assert msg.startswith("downloaded from farm.cse.ucdavis.edu")


class _FakeResponse:
    def __init__(self, chunks: list[bytes], headers: dict[str, str] | None = None):
        self._chunks = list(chunks)
        self.headers = headers or {}

    def read(self, _size: int = -1) -> bytes:
        return self._chunks.pop(0) if self._chunks else b""

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class TestDownloadFile:
    def test_streams_chunks_and_reports_progress(self, tmp_path, monkeypatch, capsys):
        response = _FakeResponse([b"abc", b"def", b""], {"Content-Length": "6"})
        monkeypatch.setattr(db.urllib.request, "urlopen", lambda *a, **k: response)
        target = tmp_path / "out.bin"
        db._download_file("https://example.invalid/f", target)
        assert target.read_bytes() == b"abcdef"
        assert "100%" in capsys.readouterr().out

    def test_unknown_length_skips_percentage(self, tmp_path, monkeypatch, capsys):
        response = _FakeResponse([b"abc", b""])
        monkeypatch.setattr(db.urllib.request, "urlopen", lambda *a, **k: response)
        db._download_file("https://example.invalid/f", tmp_path / "out.bin")
        assert "%" not in capsys.readouterr().out


class TestPanelArchiveFailure:
    def test_corrupt_panel_archive_reports(self, root, monkeypatch, tmp_path):
        bad = tmp_path / "panel.tar.gz"
        bad.write_bytes(b"not a tarball")
        (root / "D2_ani").mkdir(parents=True, exist_ok=True)
        monkeypatch.setattr(db, "_download_file", lambda url, dst: shutil.copy2(bad, dst))
        monkeypatch.setattr(db, "PANEL_SHA256", hashlib.sha256(b"not a tarball").hexdigest())
        msg = db._try_download_panel_release(root / "D2_ani")
        assert msg is not None and msg.startswith("ERROR: unsafe or corrupt archive")


class TestPanelGenomeUnzip:
    def test_zip_is_extracted_and_cleaned(self, tmp_path, monkeypatch):
        import types as _types

        dst = tmp_path / "panel"
        genomes = dst / "genomes"
        genomes.mkdir(parents=True)

        def fake_run(cmd, **kwargs):
            if cmd[0] == "/bin/datasets":
                Path(cmd[cmd.index("--filename") + 1]).write_bytes(b"zip")
            return _types.SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(db.subprocess, "run", fake_run)
        monkeypatch.setattr(
            db.subprocess,
            "run",
            lambda cmd, **kw: (
                (
                    (dst / "_dl").mkdir(exist_ok=True),
                    ((dst / "_dl" / "g_genomic.fna").write_text(">x\nACGT\n")),
                    _types.SimpleNamespace(returncode=0, stdout="", stderr=""),
                )[-1]
                if cmd[0] == "unzip"
                else fake_run(cmd, **kw)
            ),
        )
        error = db._download_panel_genomes("/bin/datasets", ["GCF_1"], dst, genomes)
        assert error is None
        assert (genomes / "g_genomic.fna").is_file()
        assert not (dst / "_genomes.zip").exists()


class TestStatusCsv:
    def test_csv_output(self, root):
        _touch(root / "D1_marker" / "markers.fasta")
        text = db.format_status("csv").splitlines()
        assert text[0] == "COMPONENT,READY,TIER,INPUTS,PATH"
        assert text[1].startswith("D1_marker,ready,mini,FASTA/FASTQ,")


class TestNamedDatabasesListing:
    def _named_marker_db(self, root):
        path = Path(root) / "marker" / "food"
        path.mkdir(parents=True)
        (path / "markers.fasta").write_text(">markers~~~g1~~~x marker\nACGTACGT\n")
        (path / "marker_rules.yaml").write_text("rules:\n- species: X\n  genes: [g1]\n")
        return path

    def test_json_includes_custom_databases(self, root):
        import json as json_mod

        self._named_marker_db(root)
        payload = json_mod.loads(db.format_tiers("json"))
        custom = payload["custom_databases"]["marker"]
        assert [entry["name"] for entry in custom] == ["food"]
        assert custom[0]["path"].endswith("marker/food")

    def test_table_appends_named_section(self, root):
        self._named_marker_db(root)
        text = db.format_tiers("table")
        assert "named databases" in text
        assert "food" in text

    def test_no_named_databases_omit_section(self, root):
        assert "named databases" not in db.format_tiers("table")
        payload = json.loads(db.format_tiers("json"))
        assert payload["custom_databases"] == {}

    def test_tsv_unchanged(self, root):
        self._named_marker_db(root)
        lines = db.format_tiers("tsv").splitlines()
        assert lines[0] == "NAME\tREADY\tSIZE\tINPUTS\tDESCRIPTION"
        assert all(len(line.split("\t")) == 5 for line in lines)


class TestFlattenIfNested:
    def test_moves_files_up_from_single_nested_dir(self, tmp_path):
        nested = tmp_path / "inner"
        (nested / "panel.sketch").mkdir(parents=True)
        (nested / "panel.sketch" / "sketches.db").write_bytes(b"x")
        db._flatten_if_nested(tmp_path, "panel.sketch/sketches.db")
        assert (tmp_path / "panel.sketch" / "sketches.db").is_file()
        assert not nested.exists()

    def test_noop_when_marker_present(self, tmp_path):
        (tmp_path / "panel.sketch").mkdir()
        (tmp_path / "panel.sketch" / "sketches.db").write_bytes(b"x")
        db._flatten_if_nested(tmp_path, "panel.sketch/sketches.db")
        assert (tmp_path / "panel.sketch" / "sketches.db").is_file()


class TestBuildFromManifestErrors:
    def test_genome_download_error_short_circuits(self, monkeypatch, tmp_path):
        manifests = tmp_path / "manifests"
        manifests.mkdir()
        (manifests / "panel_accessions.tsv").write_text("accession\torganism\nGCF_1\tX\n")
        (manifests / "metadata.tsv").write_text("fna\tspecies\n")
        monkeypatch.setattr(db, "_MANIFEST_DIR", manifests)
        monkeypatch.setattr(db, "which", lambda tool: f"/bin/{tool}")
        monkeypatch.setattr(
            db, "_download_panel_genomes", lambda *a, **k: "ERROR: NCBI datasets download failed"
        )
        assert db._build_panel_from_manifest(tmp_path / "dst").startswith("ERROR: NCBI")


class TestDownloadPanelGenomes:
    def test_ncbi_failure_is_actionable(self, monkeypatch, tmp_path):
        import types as _types

        dst = tmp_path / "dst"
        genomes = tmp_path / "genomes"
        dst.mkdir()
        genomes.mkdir()
        failed = _types.SimpleNamespace(returncode=1, stderr="quota exceeded")

        def fake_run(cmd, **kwargs):
            return failed

        monkeypatch.setattr(db.subprocess, "run", fake_run)
        error = db._download_panel_genomes("/bin/datasets", ["GCF_1"], dst, genomes)
        assert error is not None
        assert error.startswith("ERROR: NCBI datasets download failed")
