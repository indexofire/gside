"""Manifest build integration (real skani, download skipped, tiny fixture)."""

from __future__ import annotations

import shutil
import types
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
_ECOLI = _REPO / "data" / "db" / "D2_ani" / "genomes" / "GCF_000005845.2_ASM584v2_genomic.fna"

_has_skani = shutil.which("skani") is not None


@pytest.mark.skipif(not _has_skani, reason="skani not available")
@pytest.mark.skipif(not _ECOLI.is_file(), reason="panel genome missing")
class TestBuildFromManifest:
    def test_skip_download_builds_sketch(self, tmp_path, monkeypatch):
        import gside.db as dbmod
        from gside.db import _build_panel_from_manifest

        manifest = tmp_path / "manifests"
        manifest.mkdir()
        (manifest / "panel_accessions.tsv").write_text(
            "accession\torganism\nGCF_000005845.2\tEscherichia coli\n"
        )
        (manifest / "metadata.tsv").write_text("genome\tspecies\nA\tB\n")
        genomes = tmp_path / "dst" / "genomes"
        genomes.mkdir(parents=True)
        shutil.copy(_ECOLI, genomes / _ECOLI.name)
        monkeypatch.setattr(dbmod, "_MANIFEST_DIR", manifest)
        # NOTE: real skani resolves via PATH; datasets is stubbed because the
        # fixture genome is already present, so no download is attempted.
        import shutil as _shutil

        monkeypatch.setattr(
            dbmod,
            "which",
            lambda tool: _shutil.which(tool) if tool == "skani" else "/bin/datasets",
        )
        msg = _build_panel_from_manifest(tmp_path / "dst")
        assert msg.startswith("built from manifest: 1 genomes")
        assert (tmp_path / "dst" / "panel.sketch").exists()
        assert (tmp_path / "dst" / "metadata.tsv").is_file()


class TestBuildErrors:
    def test_no_genomes_after_download(self, tmp_path, monkeypatch):
        import subprocess

        import gside.db as dbmod

        manifest = tmp_path / "manifests"
        manifest.mkdir()
        (manifest / "panel_accessions.tsv").write_text("accession\nGCF_000005845.2\n")
        (manifest / "metadata.tsv").write_text("x\n")
        (tmp_path / "dst" / "genomes").mkdir(parents=True)
        monkeypatch.setattr(dbmod, "_MANIFEST_DIR", manifest)
        monkeypatch.setattr(dbmod, "which", lambda tool: "/bin/x")
        monkeypatch.setattr(
            subprocess,
            "run",
            lambda *a, **k: types.SimpleNamespace(returncode=0, stdout="", stderr=""),
        )
        assert dbmod._build_panel_from_manifest(tmp_path / "dst").startswith(
            "ERROR: no genome files"
        )

    def test_skani_failure(self, tmp_path, monkeypatch):
        import shutil
        import subprocess

        import gside.db as dbmod

        gsrc = (
            Path(__file__).resolve().parents[2]
            / "data"
            / "db"
            / "D2_ani"
            / "genomes"
            / "GCF_000005845.2_ASM584v2_genomic.fna"
        )
        if not gsrc.is_file():
            pytest.skip("genome missing")
        manifest = tmp_path / "manifests"
        manifest.mkdir()
        (manifest / "panel_accessions.tsv").write_text("accession\nGCF_000005845.2\n")
        (manifest / "metadata.tsv").write_text("x\n")
        genomes = tmp_path / "dst" / "genomes"
        genomes.mkdir(parents=True)
        shutil.copy(gsrc, genomes / gsrc.name)
        (tmp_path / "dst" / "panel.sketch").mkdir()
        monkeypatch.setattr(dbmod, "_MANIFEST_DIR", manifest)
        monkeypatch.setattr(dbmod, "which", lambda tool: "/bin/x")

        def fake_run(cmd, **kw):
            if "sketch" in cmd:
                return types.SimpleNamespace(returncode=1, stdout="", stderr="boom")
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(subprocess, "run", fake_run)
        assert dbmod._build_panel_from_manifest(tmp_path / "dst").startswith(
            "ERROR: skani sketch failed"
        )
