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

import csv
import gzip
import hashlib
import io
import json
import os
import shutil
import subprocess
import tarfile
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import Any

from rich import box
from rich.console import Console
from rich.table import Table

from .config import SPECIES_DB_DIR, which

_MANIFEST_DIR = SPECIES_DB_DIR / "D2_ani" / "manifests"

PANEL_RELEASE_URL = (
    "https://github.com/indexofire/gside/releases/download/db-v0.1/panel.sketch.tar.gz"
)
MASH_ZENODO_URL = "https://zenodo.org/records/22664519/files/RefSeqSketches_237.msh.gz"
MASH_MD5 = "dee53b23af3ab120333f9eb1b95ae60f"
# Pinned SHA256 digests (hex) of downloaded artifacts. An empty string means
# no checksum is pinned yet; _verify_sha256 prints a warning and skips.
PANEL_SHA256 = "f9bdbeeaa744a6ed9fafef98e854d309eebc8d060df787345c2d7cdd1e36e595"
# Upstream publishes no digest for the 3.9GB signature zip; this value was
# pinned after a locally verified download (2026-10-21): complete size, zip
# integrity, 143,385 k=31 DNA signatures, CC0. Override at runtime with
# GSIDE_SOURMASH_SIG_SHA256=<hex>.
SOURMASH_SIG_SHA256 = "6f3cc74f0bdcf6f84466ba91f87d6f2bf30081d5814bdf006cafd1e9cbbeada9"
SOURMASH_LINEAGES_SHA256 = "98bceab27a50f08b2f777ca7bdfb57c88aabe5ce1fa54cfb103dd0ad51b67624"
SOURMASH_FARM_BASE = "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226"
SOURMASH_SIG_URL = f"{SOURMASH_FARM_BASE}/gtdb-rs226-reps.k31.sig.zip"
SOURMASH_LINEAGES_URL = f"{SOURMASH_FARM_BASE}/gtdb-rs226-reps.lineages.csv"
# Upstream landing pages referenced by db setup / error hints.
K2_INDEXES_URL = "https://benlangmead.github.io/aws-indexes/k2"
SYLPH_DOCS_URL = "https://sylph-docs.github.io/"

# Auto-downloadable metagenomic tiers (Ben Langmead's Kraken2 packages include
# the Bracken 50-300-mer distributions; sylph ships pre-sketched GTDB DBs).
K2_STANDARD8_URL = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz"
K2_STANDARD8_MD5_URL = (
    "https://genome-idx.s3.amazonaws.com/kraken/standard_08_GB_20260626/standard_08_GB.md5"
)
SYLPH_DB_URL = "https://storage.googleapis.com/sylph-stuff/gtdb-r226-c200-dbv2.syl2db"
SYLPH_DB_FILENAME = "gtdb-r226-c200-dbv2.syl2db"
# GCS publishes the object MD5 (x-goog-hash) for this 14GB blob.
SYLPH_DB_MD5 = "1815389774987c7ce3d713c975a8853f"

_TIERS: dict[str, dict[str, Any]] = {
    "mini": {
        "description": "Marker rules + sequences (bundled)",
        "size": "~1MB",
        "items": ["D1_marker"],
    },
    "panel": {
        "description": "Curated reference panel for ANI (skani)",
        "size": "~130MB sketch / ~1GB genomes",
        "items": ["D2_ani"],
    },
    "mash": {
        "description": "RefSeq MinHash sketch",
        "size": "~179MB",
        "items": ["D3_mash"],
    },
    "sourmash": {
        "description": "GTDB gather database (sourmash, sketch k=31)",
        "size": "~3.9GB",
        "items": ["D4_sourmash"],
    },
    "meta": {
        "description": "Kraken2 Standard-8 prebuilt index (metagenomic profiling)",
        "size": "~5.5GB",
        "items": ["D5_kraken2"],
        "help_url": K2_INDEXES_URL,
    },
    "meta_sylph": {
        "description": "sylph GTDB-r226 database (ANI-based metagenomic profiling)",
        "size": "~14GB",
        "items": ["D6_sylph"],
        "help_url": SYLPH_DOCS_URL,
    },
    "all": {
        "description": "panel + mash",
        "size": "~310MB sketches",
        "items": ["D2_ani", "D3_mash"],
    },
}

_DOWNLOAD_TIMEOUT = 3600
_CHUNK = 1024 * 1024


def _has_sylph_db(db_dir: Path) -> bool:
    """True when any sylph database lives under ``db_dir`` (either format)."""
    return any(db_dir.glob("*.syldb")) or any(db_dir.glob("*.syl2db"))


# Database components in `db status` order — key: (tier, directory under
# SPECIES_DB_DIR, readiness probes relative to the directory). Ready when
# every file of one probe group is present; ``is_file`` probes require a
# regular file, ``exists`` probes match any path type.
_COMPONENTS: dict[str, tuple[str, str, list[list[tuple[str, Callable[[Path], bool]]]]]] = {
    "D1_marker": ("mini", "D1_marker", [[("markers.fasta", Path.exists)]]),
    "D2_ani": (
        "panel",
        "D2_ani",
        [[("panel.sketch/sketches.db", Path.exists)], [("panel.sketch", Path.is_file)]],
    ),
    "D3_mash": ("mash", "D3_mash", [[("mash.msh", Path.exists)], [("payload.bin", Path.exists)]]),
    "D4_sourmash": (
        "sourmash",
        "D4_sourmash",
        [[("gtdb-reps-k31.zip", Path.exists), ("lineages.csv", Path.exists)]],
    ),
    "D5_kraken2": (
        "meta",
        "D5_kraken2",
        [
            [
                ("hash.k2d", Path.exists),
                ("taxo.k2d", Path.exists),
                ("opts.k2d", Path.exists),
                # Bracken needs the distribution matching our read length (150).
                ("database150mers.kmer_distrib", Path.exists),
            ],
        ],
    ),
    "D6_sylph": ("meta_sylph", "D6_sylph", [[("", _has_sylph_db)]]),
}


def _component_ready(component: str) -> bool:
    _, dir_name, groups = _COMPONENTS[component]
    base = SPECIES_DB_DIR / dir_name
    return any(all(pred(base / rel) for rel, pred in group) for group in groups)


def _check_panel() -> bool:
    return _component_ready("D2_ani")


def _check_mash() -> bool:
    return _component_ready("D3_mash")


def _check_sourmash() -> bool:
    return _component_ready("D4_sourmash")


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

    # Installers are resolved by name at call time on purpose: tests monkeypatch
    # the module-level downloaders, so they must not be captured in a table.
    results: dict[str, str] = {}
    if tier in ("panel", "all"):
        results["D2_ani"] = _setup_component("D2_ani", _install_panel, source)
    if tier in ("mash", "all"):
        results["D3_mash"] = _setup_component("D3_mash", _download_mash_zenodo, source)
    if tier == "sourmash":
        results["D4_sourmash"] = _setup_component("D4_sourmash", _download_sourmash_farm, source)
    if tier == "meta":
        results["D5_kraken2"] = _setup_component("D5_kraken2", _install_kraken2, source)
    if tier == "meta_sylph":
        results["D6_sylph"] = _setup_component("D6_sylph", _install_sylph, source)
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


def _file_md5(path: Path) -> str:
    """Streaming MD5 of a local file (lowercase hex)."""
    digest = hashlib.md5()
    with path.open("rb") as fh:
        while chunk := fh.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def _flatten_if_nested(dst: Path, marker: str) -> None:
    """Move files up when an archive extracted into a single subdirectory."""
    if (dst / marker).exists():
        return
    nested = [p for p in dst.iterdir() if p.is_dir() and (p / marker).exists()]
    if len(nested) == 1:
        for child in nested[0].iterdir():
            shutil.move(str(child), str(dst / child.name))
        nested[0].rmdir()


def _published_md5_manifest(url: str) -> dict[str, str]:
    """``<md5>  <filename>`` lines (Kraken2 style) → {filename: md5}."""
    with urllib.request.urlopen(url, timeout=_DOWNLOAD_TIMEOUT) as resp:
        text = resp.read().decode("utf-8", "replace")
    manifest: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 32:
            manifest[parts[1].lstrip("*")] = parts[0].lower()
    return manifest


def _verify_kraken2_files(dst: Path) -> str | None:
    """Compare every extracted file against the published MD5 manifest."""
    try:
        manifest = _published_md5_manifest(K2_STANDARD8_MD5_URL)
    except Exception as e:  # noqa: BLE001 — network fallback: warn, keep the DB
        print(f"  WARNING: could not fetch Kraken2 checksums ({e}); skipping verification")
        return None
    if not manifest:
        print("  WARNING: empty Kraken2 checksum manifest; skipping verification")
        return None
    for name, expected in manifest.items():
        target = dst / name
        if not target.is_file():
            continue
        if _file_md5(target) != expected:
            return f"ERROR: MD5 mismatch for {name}"
    return None


def _install_kraken2(dst: Path) -> str:
    """Download the Kraken2 Standard-8 archive (includes Bracken distributions)."""
    archive = dst / "k2_standard_08_GB.tar.gz"
    print("  downloading Kraken2 Standard-8 (~5.5GB, includes Bracken)...")
    try:
        _download_file(K2_STANDARD8_URL, archive)
    except Exception as e:  # noqa: BLE001 — setup boundary: report, don't crash
        archive.unlink(missing_ok=True)
        return f"ERROR: Kraken2 download failed: {e}"
    print("  extracting...")
    try:
        with tarfile.open(archive, "r:gz") as tf:
            tf.extractall(dst, filter="data")
    except tarfile.TarError as e:
        archive.unlink(missing_ok=True)
        return f"ERROR: unsafe or corrupt archive: {e}"
    archive.unlink(missing_ok=True)
    _flatten_if_nested(dst, "hash.k2d")
    error = _verify_kraken2_files(dst)
    if error:
        return error
    if not _component_ready("D5_kraken2"):
        return "ERROR: archive did not contain a complete Kraken2 + Bracken database"
    return "downloaded Kraken2 Standard-8 (index + Bracken 150-mer)"


def _install_sylph(dst: Path) -> str:
    """Download a pre-sketched sylph GTDB database (-c 200, two-stage)."""
    target = dst / SYLPH_DB_FILENAME
    print("  downloading sylph GTDB r226 (~14GB, resumable by hand if it fails)...")
    try:
        _download_file(SYLPH_DB_URL, target)
    except Exception as e:  # noqa: BLE001 — setup boundary: report, don't crash
        target.unlink(missing_ok=True)
        return f"ERROR: sylph download failed: {e}"
    if SYLPH_DB_MD5 and _file_md5(target) != SYLPH_DB_MD5:
        target.unlink(missing_ok=True)
        return f"ERROR: MD5 mismatch for {target.name}"
    if not _has_sylph_db(dst):
        return "ERROR: sylph database missing after download"
    return "downloaded GTDB r226 (c200, two-stage .syl2db)"


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
    # Upstream publishes no digest for the 3.9GB signature zip; operators can
    # pin one out-of-band without editing the package.
    sig_expected = os.environ.get("GSIDE_SOURMASH_SIG_SHA256", SOURMASH_SIG_SHA256)
    for artifact, expected in (
        (sig_zip, sig_expected),
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

        digest = _file_md5(gz_path)
        if digest != MASH_MD5:
            gz_path.unlink(missing_ok=True)
            return f"ERROR: MD5 mismatch (got {digest}, expected {MASH_MD5})"

        print("  verifying MD5 ✓, decompressing...")
        with gzip.open(gz_path, "rb") as fin, open(dst / "mash.msh", "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gz_path.unlink(missing_ok=True)
        return "downloaded from Zenodo (community sketch)"
    except Exception as e:
        return f"ERROR: Zenodo download failed: {e}"


def _manifest_accessions(manifest: Path) -> list[str]:
    """GCF accessions listed in the panel manifest (header row skipped)."""
    accessions = []
    for line in manifest.read_text().splitlines()[1:]:
        acc = line.split("\t")[0].strip()
        if acc and acc.startswith("GCF_"):
            accessions.append(acc)
    return accessions


def _download_panel_genomes(
    datasets_bin: str, accessions: list[str], dst: Path, genomes_dir: Path
) -> str | None:
    """Fetch any genomes missing from ``genomes_dir``; returns an error or None."""
    existing = {f.name.split("_ASM")[0] for f in genomes_dir.glob("*.fna")}
    to_download = [a for a in accessions if a not in existing]
    if not to_download:
        return None

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
    return None


def _build_skani_sketch(skani_bin: str, fnas: list[Path], dst: Path) -> str | None:
    """Sketch every genome into ``panel.sketch``; returns an error or None."""
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
    return None


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
    error = _download_panel_genomes(datasets_bin, _manifest_accessions(manifest), dst, genomes_dir)
    if error:
        return error

    fnas = sorted(genomes_dir.glob("*.fna"))
    if not fnas:
        return "ERROR: no genome files after download"

    error = _build_skani_sketch(skani_bin, fnas, dst)
    if error:
        return error

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


def _check_format(fmt: str) -> None:
    if fmt not in ("table", "tsv", "csv", "json"):
        raise ValueError(f"unknown format {fmt!r}; expected 'table', 'tsv', 'csv', or 'json'")


def format_status(fmt: str) -> str:
    """Render ``db status`` output as table (human), tsv, or json."""
    _check_format(fmt)
    status = db_status()
    if fmt == "json":
        return json.dumps(status, indent=2, ensure_ascii=False)
    if fmt == "tsv":
        lines = ["COMPONENT\tREADY\tTIER\tPATH"]
        lines += [
            f"{name}\t{'ready' if info['ready'] else 'missing'}\t{info['tier']}\t{info['path']}"
            for name, info in status.items()
        ]
        return "\n".join(lines)
    if fmt == "csv":
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["COMPONENT", "READY", "TIER", "PATH"])
        writer.writerows(
            [name, "ready" if info["ready"] else "missing", info["tier"], info["path"]]
            for name, info in status.items()
        )
        return out.getvalue().rstrip("\n")
    return _render_table(
        "gside database status",
        ("COMPONENT", "READY", "TIER", "PATH"),
        [
            (
                name,
                "ready" if info["ready"] else "missing",
                str(info["tier"]),
                str(info["path"]),
            )
            for name, info in status.items()
        ],
        caption="Use: gside db setup --tier <tier>",
    )


def format_tiers(fmt: str) -> str:
    """Render ``db list`` output as table (human), tsv, csv, or json."""
    _check_format(fmt)
    status = db_status()
    rows = [
        (
            name,
            "ready" if all(status[item]["ready"] for item in info["items"]) else "missing",
            info["size"],
            info["description"],
        )
        for name, info in _TIERS.items()
    ]
    if fmt == "json":
        payload = {
            name: {**info, "ready": ready == "ready"}
            for (name, ready, _, _), (_, info) in zip(rows, _TIERS.items(), strict=True)
        }
        return json.dumps(payload, indent=2, ensure_ascii=False)
    if fmt == "tsv":
        lines = ["NAME\tREADY\tSIZE\tDESCRIPTION"]
        lines += ["\t".join(row) for row in rows]
        return "\n".join(lines)
    if fmt == "csv":
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["NAME", "READY", "SIZE", "DESCRIPTION"])
        writer.writerows(rows)
        return out.getvalue().rstrip("\n")
    return _render_table(
        "gside database tiers",
        ("TIER", "READY", "SIZE", "DESCRIPTION"),
        rows,
        caption="Tiers marked ready can be used without downloading.",
    )


def _render_table(
    title: str,
    columns: tuple[str, ...],
    rows: list[tuple[str, str, str, str]],
    *,
    caption: str | None = None,
) -> str:
    table = Table(
        title=title,
        caption=caption,
        box=box.SQUARE,
        header_style="bold cyan",
        show_lines=False,
        padding=(0, 1),
    )
    for column in columns:
        table.add_column(column)
    for row in rows:
        table.add_row(*row)
    buffer = io.StringIO()
    Console(
        file=buffer,
        force_terminal=False,
        color_system=None,
        width=120,
    ).print(table)
    return buffer.getvalue().rstrip("\n")


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
        _print_tiers()
        return 0
    else:
        print(f"unknown subcommand: {sub} (status|setup|list)")
        return 1


def _print_status() -> None:
    print(format_status("table"))


def _print_tiers() -> None:
    print(format_tiers("table"))
