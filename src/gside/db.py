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
import shutil
import subprocess
import tarfile
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import Any

from .config import SPECIES_DB_DIR, which

_MANIFEST_DIR = SPECIES_DB_DIR / "L2_ani" / "manifests"

PANEL_RELEASE_URL = (
    "https://github.com/indexofire/gside/releases/download/db-v0.1/panel.sketch.tar.gz"
)
MASH_ZENODO_URL = "https://zenodo.org/records/22664519/files/RefSeqSketches_237.msh.gz"
MASH_MD5 = "dee53b23af3ab120333f9eb1b95ae60f"
# Pinned SHA256 digests (hex) of downloaded artifacts. An empty string means
# no checksum is pinned yet; _verify_sha256 prints a warning and skips.
PANEL_SHA256 = "f9bdbeeaa744a6ed9fafef98e854d309eebc8d060df787345c2d7cdd1e36e595"
SOURMASH_SIG_SHA256 = ""
SOURMASH_LINEAGES_SHA256 = "98bceab27a50f08b2f777ca7bdfb57c88aabe5ce1fa54cfb103dd0ad51b67624"
SOURMASH_FARM_BASE = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226"
SOURMASH_SIG_URL = f"{SOURMASH_FARM_BASE}/gtdb-rs226-reps.k31.sig.zip"
SOURMASH_LINEAGES_URL = f"{SOURMASH_FARM_BASE}/gtdb-rs226-reps.lineages.csv"

_TIERS: dict[str, dict[str, Any]] = {
    "mini": {
        "description": "Marker rules + sequences (bundled)",
        "size": "~1MB",
        "items": ["markers"],
    },
    "panel": {
        "description": "Curated reference panel for ANI (skani)",
        "size": "~130MB sketch / ~1GB genomes",
        "items": ["L2_ani"],
    },
    "mash": {
        "description": "RefSeq MinHash sketch",
        "size": "~179MB",
        "items": ["L3_mash"],
    },
    "sourmash": {
        "description": "GTDB gather database (sourmash, sketch k=31)",
        "size": "~3.9GB",
        "items": ["L4_sourmash"],
    },
    "all": {
        "description": "panel + mash",
        "size": "~310MB sketches",
        "items": ["L2_ani", "L3_mash"],
    },
}

_DOWNLOAD_TIMEOUT = 3600
_CHUNK = 1024 * 1024


# Database components in `db status` order — key: (tier, directory under
# SPECIES_DB_DIR, readiness probes relative to the directory). Ready when
# every file of one probe group is present; ``is_file`` probes require a
# regular file, ``exists`` probes match any path type.
_COMPONENTS: dict[str, tuple[str, str, list[list[tuple[str, Callable[[Path], bool]]]]]] = {
    "markers": ("mini", "L1_marker", [[("markers.fasta", Path.exists)]]),
    "L2_ani": (
        "panel",
        "L2_ani",
        [[("panel.sketch/sketches.db", Path.exists)], [("panel.sketch", Path.is_file)]],
    ),
    "L3_mash": ("mash", "L3_mash", [[("mash.msh", Path.exists)], [("payload.bin", Path.exists)]]),
    "L4_sourmash": (
        "sourmash",
        "L4_sourmash",
        [[("gtdb-reps-k31.zip", Path.exists), ("lineages.csv", Path.exists)]],
    ),
}


def _component_ready(component: str) -> bool:
    _, dir_name, groups = _COMPONENTS[component]
    base = SPECIES_DB_DIR / dir_name
    return any(all(pred(base / rel) for rel, pred in group) for group in groups)


def _check_panel() -> bool:
    return _component_ready("L2_ani")


def _check_mash() -> bool:
    return _component_ready("L3_mash")


def _check_sourmash() -> bool:
    return _component_ready("L4_sourmash")


def db_status() -> dict[str, dict[str, Any]]:
    return {
        component: {
            "ready": _component_ready(component),
            "path": str(SPECIES_DB_DIR / dir_name),
            "tier": tier,
        }
        for component, (tier, dir_name, _probes) in _COMPONENTS.items()
    }


def db_setup(tier: str = "panel", source: str = "") -> dict[str, str]:
    if tier not in _TIERS:
        raise ValueError(f"unknown tier {tier!r}; expected one of {sorted(_TIERS)}")

    results: dict[str, str] = {}
    if tier in ("panel", "all"):
        results["L2_ani"] = _setup_component("L2_ani", _install_panel, source)
    if tier in ("mash", "all"):
        results["L3_mash"] = _setup_component("L3_mash", _download_mash_zenodo, source)
    if tier == "sourmash":
        results["L4_sourmash"] = _setup_component("L4_sourmash", _download_sourmash_farm, source)
    if not results:
        results["info"] = "mini tier ships with repo"
    return results


def _setup_component(component: str, install: Callable[[Path], str], source: str = "") -> str:
    if _component_ready(component):
        return "already ready"

    dir_name = _COMPONENTS[component][1]
    if source:
        return _copy_from_source(source, dir_name)

    dst = SPECIES_DB_DIR / dir_name
    dst.mkdir(parents=True, exist_ok=True)
    return install(dst)


def _install_panel(dst: Path) -> str:
    """Panel install chain: pre-built GitHub Release sketch, else manifest build."""
    result = _try_download_panel_release(dst)
    return result if result else _build_panel_from_manifest(dst)


def _download_sourmash_farm(dst: Path) -> str:
    # No checksum published upstream for the sig zip; the lineages csv is pinned.
    try:
        print("  downloading sourmash GTDB sketch (~3.9GB)...")
        sig_zip = dst / "gtdb-rs226-reps.k31.sig.zip"
        lineages_csv = dst / "gtdb-rs226-reps.lineages.csv"
        _download_file(SOURMASH_SIG_URL, sig_zip)
        print("  downloading GTDB lineages...")
        _download_file(SOURMASH_LINEAGES_URL, lineages_csv)
    except Exception as e:
        return f"ERROR: farm download failed: {e}"
    for artifact, expected in (
        (sig_zip, SOURMASH_SIG_SHA256),
        (lineages_csv, SOURMASH_LINEAGES_SHA256),
    ):
        if not _verify_sha256(artifact, expected):
            artifact.unlink(missing_ok=True)
            return f"ERROR: SHA256 mismatch for {artifact.name}"
    try:
        sig_zip.rename(dst / "gtdb-reps-k31.zip")
        lineages_csv.rename(dst / "lineages.csv")
    except OSError as e:
        return f"ERROR: {e}"
    if _check_sourmash():
        return "downloaded from farm.cse.ucdavis.edu (GTDB rs226 reps, renamed in place)"
    return "ERROR: download incomplete, expected gtdb-reps-k31.zip + lineages.csv"


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
        print("  downloading panel sketch from GitHub Release...")
        archive = dst / "panel.sketch.tar.gz"
        _download_file(PANEL_RELEASE_URL, archive)
        if not _verify_sha256(archive, PANEL_SHA256):
            archive.unlink(missing_ok=True)
            return f"ERROR: SHA256 mismatch for {archive.name}"
        print("  extracting...")
        try:
            with tarfile.open(archive, "r:gz") as tf:
                tf.extractall(dst, filter="data")
        except tarfile.TarError as e:
            archive.unlink(missing_ok=True)
            return f"ERROR: unsafe or corrupt archive: {e}"
        archive.unlink(missing_ok=True)
        if _check_panel():
            return "downloaded from GitHub Release (pre-built sketch)"
        return "ERROR: downloaded archive did not contain valid panel.sketch"
    except Exception as e:
        print(f"  GitHub Release download failed: {e}")
        return None


def _download_mash_zenodo(dst: Path) -> str:
    try:
        print("  downloading mash sketch from Zenodo (~179MB)...")
        gz_path = dst / "mash.msh.gz"
        _download_file(MASH_ZENODO_URL, gz_path)

        md5 = hashlib.md5()
        with gz_path.open("rb") as f:
            while chunk := f.read(_CHUNK):
                md5.update(chunk)
        if md5.hexdigest() != MASH_MD5:
            gz_path.unlink(missing_ok=True)
            return f"ERROR: MD5 mismatch (got {md5.hexdigest()}, expected {MASH_MD5})"

        print("  verifying MD5 ✓, decompressing...")
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
            [
                datasets_bin,
                "download",
                "genome",
                "accession",
                "-f",
                str(acc_file),
                "--include",
                "genome",
                "--filename",
                str(dst / "_genomes.zip"),
            ],
            capture_output=True,
            text=True,
            timeout=_DOWNLOAD_TIMEOUT,
        )
        acc_file.unlink(missing_ok=True)
        if result.returncode != 0:
            return f"ERROR: NCBI datasets download failed: {result.stderr[:200]}"

        zip_path = dst / "_genomes.zip"
        if zip_path.exists():
            subprocess.run(
                ["unzip", "-oq", str(zip_path), "-d", str(dst / "_dl")],
                capture_output=True,
                timeout=300,
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
        capture_output=True,
        text=True,
        timeout=_DOWNLOAD_TIMEOUT,
    )
    if result.returncode != 0:
        return f"ERROR: skani sketch failed: {result.stderr[:200]}"

    shutil.copy2(metadata, dst / "metadata.tsv")

    return f"built from manifest: {len(fnas)} genomes → panel.sketch"


def _download_file(url: str, dst: Path) -> None:
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


def _verify_sha256(path: Path, expected: str) -> bool:
    """Streaming SHA256 check; empty ``expected`` warns and skips verification."""
    if not expected:
        print(f"  WARNING: no pinned checksum for {path.name}; skipping integrity verification")
        return True
    digest = hashlib.sha256()
    try:
        with path.open("rb") as f:
            while chunk := f.read(_CHUNK):
                digest.update(chunk)
    except OSError:
        return False
    return digest.hexdigest() == expected.lower()


def run_db_command(args: list[str]) -> int:
    if not args or args[0] == "status":
        _print_status()
        return 0

    sub = args[0]
    if sub == "setup":
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
