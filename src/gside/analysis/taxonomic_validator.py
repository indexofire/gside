"""Assembly quality + taxonomy validation (NCBI-style: CheckM2 + GTDB-Tk).

`gside validate` answers a different question from `gside scan`: "is this
assembly good enough, and what does GTDB-Tk place it as?" — not "which
species is this?". Both tools are required for a complete report:

* CheckM2 → completeness / contamination (quality verdict)
* GTDB-Tk → GTDB classification

Missing databases or binaries are recorded per tool (``status`` ok/partial/
error) instead of crashing, and every run writes ``validation.json`` with the
genome name, tool versions, parameters and results.
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from gside import __version__
from gside.config import CHECKM2_DB, GTDB_DB, pixi_path

logger = logging.getLogger(__name__)

SCHEMA = "gside.validate/1"
_CHECKM2_TIMEOUT = 600
_GTDBTK_TIMEOUT = 3600
_VERSION_TIMEOUT = 30


@dataclass
class ValidateResult:
    schema: str = SCHEMA
    genome: dict[str, str] = field(default_factory=dict)
    parameters: dict[str, Any] = field(default_factory=dict)
    software: dict[str, str] = field(default_factory=dict)
    tools: dict[str, str] = field(default_factory=dict)
    quality: dict[str, Any] = field(default_factory=dict)
    taxonomy: dict[str, str] = field(default_factory=dict)
    status: str = "ok"
    errors: list[str] = field(default_factory=list)
    interpretation: str = ""
    created_at: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "genome": self.genome,
            "parameters": self.parameters,
            "software": self.software,
            "tools": self.tools,
            "quality": self.quality,
            "taxonomy": self.taxonomy,
            "status": self.status,
            "errors": self.errors,
            "interpretation": self.interpretation,
            "created_at": self.created_at,
        }


def _find_tool(name: str) -> str | None:
    return shutil.which(name, path=pixi_path())


def _tool_version(binary: str) -> str:
    """Best-effort ``<tool> --version``; empty string when unavailable."""
    try:
        proc = subprocess.run(
            [binary, "--version"],
            capture_output=True,
            text=True,
            timeout=_VERSION_TIMEOUT,
            env={**__import__("os").environ, "PATH": pixi_path()},
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = (proc.stdout or proc.stderr or "").strip()
    return text.splitlines()[0].strip() if text else ""


def _parse_tsv_row(path: Path) -> dict[str, str]:
    lines = [line for line in path.read_text().strip().split("\n") if line.strip()]
    if len(lines) < 2:
        return {}
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    return dict(zip(header, values, strict=False))


def _as_float(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _run_checkm2(contigs_path: str, output_dir: Path, threads: int = 4) -> dict[str, Any]:
    """Run CheckM2; returns status/completeness/contamination/version/error."""
    result: dict[str, Any] = {
        "status": "missing",
        "completeness": None,
        "contamination": None,
        "version": "",
        "error": "",
    }
    if not CHECKM2_DB:
        result["error"] = "CheckM2 database not configured (set CHECKM2DB)"
        return result
    checkm2 = _find_tool("checkm2")
    if not checkm2:
        result["error"] = "checkm2 not found on PATH"
        return result

    result["version"] = _tool_version(checkm2)
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / "checkm2_results.tsv"
    cmd = [
        checkm2,
        "predict",
        "--input",
        contigs_path,
        "--output",
        str(tsv_path),
        "--database_path",
        str(CHECKM2_DB),
        "--threads",
        str(threads),
    ]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=_CHECKM2_TIMEOUT,
            env={**__import__("os").environ, "PATH": pixi_path()},
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        result.update(status="failed", error=f"CheckM2 failed: {exc}")
        return result
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()[:300]
        result.update(status="failed", error=f"CheckM2 failed (exit {proc.returncode}): {detail}")
        return result
    if not tsv_path.is_file():
        result.update(status="failed", error="CheckM2 produced no checkm2_results.tsv")
        return result

    row = _parse_tsv_row(tsv_path)
    result["completeness"] = _as_float(row.get("Completeness"))
    result["contamination"] = _as_float(row.get("Contamination"))
    if result["completeness"] is None or result["contamination"] is None:
        result.update(
            status="failed",
            error="CheckM2 report has no numeric Completeness/Contamination",
        )
        return result
    result["status"] = "ok"
    return result


def _run_gtdbtk(contigs_path: str, output_dir: Path, threads: int = 4) -> dict[str, Any]:
    """Run GTDB-Tk classify_wf; returns status/classification/version/note."""
    result: dict[str, Any] = {
        "status": "missing",
        "classification": "",
        "version": "",
        "note": "",
    }
    if not GTDB_DB:
        result["note"] = "GTDB-Tk database not configured (set GTDBTK_DATA_PATH or GTDBDB)"
        return result
    gtdbtk = _find_tool("gtdbtk")
    if not gtdbtk:
        result["note"] = "gtdbtk not found on PATH"
        return result

    result["version"] = _tool_version(gtdbtk)
    genome_dir = output_dir / "gtdb_genomes"
    genome_dir.mkdir(parents=True, exist_ok=True)
    genome_file = genome_dir / "query.fasta"
    genome_file.write_text(Path(contigs_path).read_text())
    out_dir = output_dir / "gtdb_output"
    cmd = [
        gtdbtk,
        "classify_wf",
        "--genome_dir",
        str(genome_dir),
        "--out_dir",
        str(out_dir),
        "--database",
        str(GTDB_DB),
        "--cpus",
        str(threads),
    ]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=_GTDBTK_TIMEOUT,
            env={**__import__("os").environ, "PATH": pixi_path()},
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        result.update(status="failed", note=f"GTDB-Tk failed: {exc}")
        return result
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()[:300]
        result.update(status="failed", note=f"GTDB-Tk failed (exit {proc.returncode}): {detail}")
        return result

    summary = out_dir / "gtdbtk.bac120.summary.tsv"
    if not summary.is_file():
        summary = out_dir / "gtdbtk.ar122.summary.tsv"
    if not summary.is_file():
        result.update(status="failed", note="GTDB-Tk produced no bac120/ar122 summary")
        return result
    row = _parse_tsv_row(summary)
    classification = row.get("classification", "").strip()
    if not classification:
        result.update(status="failed", note="GTDB-Tk summary has no classification")
        return result
    result.update(status="ok", classification=classification)
    return result


def _quality_verdict(completeness: float | None, contamination: float | None) -> str:
    if completeness is None or contamination is None:
        return "UNKNOWN"
    if completeness >= 90 and contamination <= 5:
        return "PASS"
    if completeness < 50 or contamination > 10:
        return "WARNING"
    return "CHECK"


def _build_interpretation(result: ValidateResult) -> str:
    parts: list[str] = []
    verdict = result.quality.get("verdict", "UNKNOWN")
    completeness = result.quality.get("completeness")
    contamination = result.quality.get("contamination")
    if completeness is not None and contamination is not None:
        parts.append(
            f"Genome quality: {verdict} (completeness={completeness:.1f}%,"
            f" contamination={contamination:.1f}%)"
        )
    else:
        parts.append(f"Genome quality: {verdict} (CheckM2 unavailable)")

    classification = result.taxonomy.get("classification", "")
    if classification:
        parts.append(f"GTDB-Tk: {classification}")
    else:
        parts.append("GTDB-Tk: unavailable")
    if result.status != "ok":
        parts.append(f"Status: {result.status}")
    return "; ".join(parts)


def validate_genome(
    contigs_path: str | Path,
    output_dir: str | Path,
    threads: int = 4,
) -> ValidateResult:
    """Run CheckM2 + GTDB-Tk and write ``<output_dir>/validation.json``."""
    contigs = Path(contigs_path).resolve()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    checkm2 = _run_checkm2(str(contigs), out, threads)
    gtdbtk = _run_gtdbtk(str(contigs), out, threads)

    quality = {
        "completeness": checkm2["completeness"],
        "contamination": checkm2["contamination"],
        "verdict": _quality_verdict(checkm2["completeness"], checkm2["contamination"]),
    }
    errors = [msg for msg in (checkm2["error"], gtdbtk["note"]) if msg]
    ok_tools = sum(1 for tool in (checkm2, gtdbtk) if tool["status"] == "ok")
    status = "ok" if ok_tools == 2 else ("partial" if ok_tools == 1 else "error")

    result = ValidateResult(
        genome={"name": contigs.stem, "path": str(contigs)},
        parameters={
            "threads": threads,
            "checkm2_db": str(CHECKM2_DB) if CHECKM2_DB else None,
            "gtdbtk_db": str(GTDB_DB) if GTDB_DB else None,
        },
        software={
            "gside": __version__,
            "checkm2": checkm2["version"],
            "gtdbtk": gtdbtk["version"],
        },
        tools={"checkm2": checkm2["status"], "gtdbtk": gtdbtk["status"]},
        quality=quality,
        taxonomy={"classification": gtdbtk["classification"], "note": gtdbtk["note"]},
        status=status,
        errors=errors,
        created_at=datetime.now(UTC).isoformat(timespec="seconds"),
    )
    result.interpretation = _build_interpretation(result)
    (out / "validation.json").write_text(
        json.dumps(result.to_dict(), ensure_ascii=False, indent=2) + "\n"
    )
    return result
