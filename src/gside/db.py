"""Database management — download, build, and status for gside reference DBs.

Download paths (in priority order):
  1. GitHub Release: pre-built sketch (fast, ~130MB)
  2. Zenodo: community mash sketch (~179MB)
  3. Build from manifest: download genomes from NCBI + build skani sketch

Tiers:
  mini     marker rules + sequences (ships with repo)
  panel    curated reference panel (skani sketch)
  mash     RefSeq MinHash sketch
  all      panel + mash
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Any

from .config import DATA_DIR, SPECIES_DB_DIR, which

_MANIFEST_DIR = DATA_DIR / "panel_manifest"

PANEL_RELEASE_URL = (
    "https://github.com/indexofire/gside/releases/download/db-v0.1/panel.sketch.tar.gz"
)
MASH_ZENODO_URL = "https://zenodo.org/records/22664519/files/RefSeqSketches_237.msh.gz"
MASH_MD5 = "dee53b23af3ab120333f9eb1b95ae60f"

_TIERS: dict[str, dict[str, Any]] = {
    "mini": {
        "description": "Marker rules + sequences (bundled)",
        "size": "~1MB",
        "items": ["markers_v2"],
    },
    "panel": {
        "description": "Curated reference panel for ANI (skani)",
        "size": "~130MB sketch / ~1GB genomes",
        "items": ["refseq_panel"],
    },
    "mash": {
        "description": "RefSeq MinHash sketch",
        "size": "~179MB",
        "items": ["mash_refseq"],
    },
    "all": {
        "description": "panel + mash",
        "size": "~310MB sketches",
        "items": ["refseq_panel", "mash_refseq"],
    },
}

_DOWNLOAD_TIMEOUT = 3600
_CHUNK = 1024 * 1024


def db_status() -> dict[str, dict[str, Any]]:
    return {
        "markers_v2": {
            "ready": (DATA_DIR / "reference/species/markers_v2.fasta").exists(),
            "path": str(DATA_DIR / "reference/species"),
            "tier": "mini",
        },
        "refseq_panel": {
            "ready": _check_panel(),
            "path": str(SPECIES_DB_DIR / "refseq_panel"),
            "tier": "panel",
        },
        "mash_refseq": {
            "ready": _check_mash(),
            "path": str(SPECIES_DB_DIR / "mash_refseq"),
            "tier": "mash",
        },
    }


def _check_panel() -> bool:
    p = SPECIES_DB_DIR / "refseq_panel"
    return (p / "panel.sketch" / "sketches.db").exists() or (
        p / "panel.sketch"
    ).is_file()


def _check_mash() -> bool:
    m = SPECIES_DB_DIR / "mash_refseq"
    return (m / "mash.msh").exists() or (m / "payload.bin").exists()


def db_setup(tier: str = "panel", source: str = "") -> dict[str, str]:
    if tier not in _TIERS:
        raise ValueError(f"unknown tier {tier!r}; expected one of {sorted(_TIERS)}")

    results: dict[str, str] = {}
    if tier in ("panel", "all"):
        results["refseq_panel"] = _setup_panel(source)
    if tier in ("mash", "all"):
        results["mash_refseq"] = _setup_mash(source)
    if not results:
        results["info"] = "mini tier ships with repo"
    return results


def _setup_panel(source: str = "") -> str:
    if _check_panel():
        return "already ready"

    if source:
        return _copy_from_source(source, "refseq_panel")

    dst = SPECIES_DB_DIR / "refseq_panel"
    dst.mkdir(parents=True, exist_ok=True)

    result = _try_download_panel_release(dst)
    if result:
        return result

    return _build_panel_from_manifest(dst)


def _setup_mash(source: str = "") -> str:
    if _check_mash():
        return "already ready"

    if source:
        return _copy_from_source(source, "mash_refseq")

    dst = SPECIES_DB_DIR / "mash_refseq"
    dst.mkdir(parents=True, exist_ok=True)

    return _download_mash_zenodo(dst)


def _copy_from_source(source: str, name: str) -> str:
    src = Path(source)
    if src.name != name:
        src = src / name
    if not src.exists():
        return f"ERROR: source not found: {src}"
    dst = SPECIES_DB_DIR / name
    dst.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst, dirs_exist_ok=True)
    return f"copied from {src}"


def _try_download_panel_release(dst: Path) -> str | None:
    try:
        print(f"  downloading panel sketch from GitHub Release...")
        archive = dst / "panel.sketch.tar.gz"
        _download_file(PANEL_RELEASE_URL, archive)
        print(f"  extracting...")
        with tarfile.open(archive, "r:gz") as tf:
            tf.extractall(dst)
        archive.unlink(missing_ok=True)
        if _check_panel():
            return "downloaded from GitHub Release (pre-built sketch)"
        return "ERROR: downloaded archive did not contain valid panel.sketch"
    except Exception as e:
        print(f"  GitHub Release download failed: {e}")
        return None


def _download_mash_zenodo(dst: Path) -> str:
    try:
        print(f"  downloading mash sketch from Zenodo (~179MB)...")
        gz_path = dst / "mash.msh.gz"
        _download_file(MASH_ZENODO_URL, gz_path)

        md5 = hashlib.md5(gz_path.read_bytes()).hexdigest()
        if md5 != MASH_MD5:
            gz_path.unlink(missing_ok=True)
            return f"ERROR: MD5 mismatch (got {md5}, expected {MASH_MD5})"

        print(f"  verifying MD5 ✓, decompressing...")
        with gzip.open(gz_path, "rb") as fin, open(dst / "mash.msh", "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gz_path.unlink(missing_ok=True)
        return "downloaded from Zenodo (community sketch)"
    except Exception as e:
        return f"ERROR: Zenodo download failed: {e}"


def _build_panel_from_manifest(dst: Path) -> str:
    manifest = _MANIFEST_DIR / "panel_accessions.tsv"
    metadata = _MANIFEST_DIR / "metadata.tsv"

    if not manifest.exists() or not metadata.exists():
        return "ERROR: panel manifest not found — cannot build from source"

    datasets_bin = which("datasets")
    if not datasets_bin:
        return "ERROR: NCBI datasets CLI not found (pixi add ncbi-datasets-cli)"

    skani_bin = which("skani")
    if not skani_bin:
        return "ERROR: skani not found (pixi add skani)"

    genomes_dir = dst / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)

    accessions = []
    for line in manifest.read_text().splitlines()[1:]:
        acc = line.split("\t")[0].strip()
        if acc and acc.startswith("GCF_"):
            accessions.append(acc)

    existing = {f.name.split("_ASM")[0] for f in genomes_dir.glob("*.fna")}
    to_download = [a for a in accessions if a not in existing]

    if to_download:
        print(f"  downloading {len(to_download)} genomes from NCBI...")
        acc_file = dst / "_accessions.txt"
        acc_file.write_text("\n".join(to_download))
        result = subprocess.run(
            [datasets_bin, "download", "genome", "accession", "-f", str(acc_file),
             "--include", "genome", "--filename", str(dst / "_genomes.zip")],
            capture_output=True, text=True, timeout=_DOWNLOAD_TIMEOUT,
        )
        acc_file.unlink(missing_ok=True)
        if result.returncode != 0:
            return f"ERROR: NCBI datasets download failed: {result.stderr[:200]}"

        zip_path = dst / "_genomes.zip"
        if zip_path.exists():
            subprocess.run(
                ["unzip", "-oq", str(zip_path), "-d", str(dst / "_dl")],
                capture_output=True, timeout=300,
            )
            for fna in (dst / "_dl").rglob("*_genomic.fna"):
                shutil.copy2(fna, genomes_dir / fna.name)
            shutil.rmtree(dst / "_dl", ignore_errors=True)
            zip_path.unlink(missing_ok=True)

    fnas = sorted(genomes_dir.glob("*.fna"))
    if not fnas:
        return "ERROR: no genome files after download"

    print(f"  building skani sketch from {len(fnas)} genomes...")
    sketch_dir = dst / "panel.sketch"
    if sketch_dir.exists():
        shutil.rmtree(sketch_dir)
    result = subprocess.run(
        [skani_bin, "sketch", *[str(f) for f in fnas], "-o", str(sketch_dir)],
        capture_output=True, text=True, timeout=_DOWNLOAD_TIMEOUT,
    )
    if result.returncode != 0:
        return f"ERROR: skani sketch failed: {result.stderr[:200]}"

    shutil.copy2(metadata, dst / "metadata.tsv")

    return f"built from manifest: {len(fnas)} genomes → panel.sketch"


def _download_file(url: str, dst: Path) -> None:
    import urllib.error

    req = urllib.request.Request(url, headers={"User-Agent": "gside-setup"})
    with urllib.request.urlopen(req, timeout=_DOWNLOAD_TIMEOUT) as resp:
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        with open(dst, "wb") as f:
            while True:
                chunk = resp.read(_CHUNK)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if total > 0:
                    pct = downloaded * 100 // total
                    print(f"\r  {pct:3d}% ({downloaded // _CHUNK}MB)", end="", flush=True)
        print()


def run_db_command(args: list[str]) -> int:
    if not args or args[0] == "status":
        _print_status()
        return 0

    sub = args[0]
    if sub == "status":
        _print_status()
        return 0
    elif sub == "setup":
        tier = "panel"
        source = ""
        i = 1
        while i < len(args):
            if args[i] == "--tier" and i + 1 < len(args):
                tier = args[i + 1]
                i += 2
            elif args[i] == "--source" and i + 1 < len(args):
                source = args[i + 1]
                i += 2
            else:
                i += 1
        results = db_setup(tier, source)
        has_error = False
        for name, msg in results.items():
            print(f"  {name}: {msg}")
            if "ERROR" in msg:
                has_error = True
        return 1 if has_error else 0
    elif sub == "list":
        for name, info in _TIERS.items():
            print(f"  {name:8} {info['size']:28} {info['description']}")
        return 0
    else:
        print(f"unknown subcommand: {sub} (status|setup|list)")
        return 1


def _print_status() -> None:
    status = db_status()
    print("gside database status")
    print("─" * 55)
    for name, info in status.items():
        icon = "✅" if info["ready"] else "❌"
        print(f"  {icon} {name:16} tier={info['tier']:6} path={info['path']}")
    print()
    print("  Run 'gside db setup --tier <tier>' to provision")
