"""Command-line interface for gside."""

from __future__ import annotations

import contextlib
import json
import logging
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any, TextIO

import click

from gside import __version__

HELP_SETTINGS = {"help_option_names": ["-h", "--help"]}

MODES = ["marker", "panel", "mash_refseq", "sourmash", "gtdbtk", "all"]
TIERS = ["mini", "panel", "mash", "sourmash", "all"]

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Configure root logging. Call once from the CLI entry point."""
    if verbose and quiet:
        raise ValueError("verbose and quiet cannot be enabled together")
    level = logging.DEBUG if verbose else logging.ERROR if quiet else logging.WARNING
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


@click.group(
    context_settings=HELP_SETTINGS,
    invoke_without_command=True,
    no_args_is_help=True,
)
@click.option("-V", "--version", is_flag=True, help="Show the version and exit.")
@click.option(
    "-J", "--json", "version_json", is_flag=True, help="With --version: emit version JSON."
)
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
@click.option("--quiet", "-q", is_flag=True, help="Suppress non-error logging.")
@click.pass_context
def main(ctx: click.Context, version: bool, version_json: bool, verbose: bool, quiet: bool) -> None:
    """gside — Genome Species IDentification Engine.

    Multi-layer species identification from contigs. Prints a single
    JSON verdict per run to stdout.
    """
    if version:
        if version_json:
            print(json.dumps({"tool": "gside", "version": __version__}))
        else:
            print(f"gside, version {__version__}")
        sys.exit(0)
    if verbose and quiet:
        raise click.UsageError("--verbose and --quiet cannot be used together")
    setup_logging(verbose=verbose, quiet=quiet)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@main.command("species", context_settings=HELP_SETTINGS)
@click.argument("contigs", nargs=-1, required=True, type=click.Path(dir_okay=False))
@click.option(
    "-m",
    "--mode",
    type=click.Choice(MODES),
    default="marker",
    show_default=True,
    help="Identification method (all runs every method, then arbitrates).",
)
@click.option(
    "-D",
    "--db-dir",
    type=click.Path(file_okay=False),
    default=None,
    help="Database root (default: $GSIDE_DB_DIR or data/db).",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["json", "tsv", "md"]),
    default="json",
    show_default=True,
    help="Output format (tsv/md render one row per input file).",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help=(
        "Write the report to PATH instead of stdout (truncates; tsv/md rows"
        " stream to the file as each input completes; json is written once"
        " at the end; stdout then carries no data)"
    ),
)
@click.option(
    "-j",
    "--jobs",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Number of input files to identify concurrently (output order is always input order).",
)
@click.option(
    "-t",
    "--threads",
    type=click.IntRange(min=1),
    default=4,
    show_default=True,
    help="Threads per aligner invocation (blastn -num_threads / skani -t / mash -p).",
)
@click.option(
    "-n",
    "--no-header",
    is_flag=True,
    default=False,
    help="Suppress the TSV header row (tsv only; md/json keep their structure).",
)
@click.option(
    "-a",
    "--auto-setup",
    is_flag=True,
    default=False,
    help="Automatically install missing database tiers before identification (progress on stderr).",
)
def species_cmd(
    contigs: tuple[str, ...],
    mode: str,
    db_dir: str | None,
    output_format: str,
    output: str | None,
    jobs: int,
    threads: int,
    no_header: bool,
    auto_setup: bool,
) -> None:
    """Identify species from assembled contigs FASTA (accepts multiple files)."""
    sys.exit(
        _run_species(
            list(contigs), mode, db_dir, output_format, output, jobs, threads, no_header, auto_setup
        )
    )


@click.group(
    "db",
    context_settings=HELP_SETTINGS,
    invoke_without_command=True,
)
@click.pass_context
def db_group(ctx: click.Context) -> None:
    """Manage reference databases (status/setup/list)."""
    if ctx.invoked_subcommand is None:
        from .db import run_db_command

        sys.exit(run_db_command(["status"]))


@db_group.command("status", context_settings=HELP_SETTINGS)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["table", "tsv", "csv", "json"]),
    default="table",
    show_default=True,
    help="Output format.",
)
def db_status_cmd(output_format: str) -> None:
    """Show database readiness status."""
    from .db import format_status

    print(format_status(output_format))


@db_group.command("setup", context_settings=HELP_SETTINGS)
@click.option(
    "-t",
    "--tier",
    type=click.Choice(TIERS),
    default="panel",
    show_default=True,
    help="Database tier to install (all = panel + mash; sourmash is separate).",
)
@click.option(
    "-s",
    "--source",
    default="",
    help="Copy from an existing local database directory.",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
    help="Output format (json routes progress to stderr, result to stdout).",
)
def db_setup_cmd(tier: str, source: str, output_format: str) -> None:
    """Install or update databases."""
    from .db import db_setup, run_db_command

    if output_format == "json":
        with contextlib.redirect_stdout(sys.stderr):
            results = db_setup(tier, source)
        print(json.dumps(results, indent=2, ensure_ascii=False))
        sys.exit(1 if any("ERROR" in str(msg) for msg in results.values()) else 0)

    args = ["setup", "--tier", tier]
    if source:
        args += ["--source", source]
    sys.exit(run_db_command(args))


@db_group.command("list", context_settings=HELP_SETTINGS)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["table", "tsv", "csv", "json"]),
    default="table",
    show_default=True,
    help="Output format.",
)
def db_list_cmd(output_format: str) -> None:
    """List available database tiers."""
    from .db import format_tiers

    print(format_tiers(output_format))


@main.command("validate", context_settings=HELP_SETTINGS)
@click.argument("contigs", type=click.Path(dir_okay=False))
@click.option(
    "--mode",
    type=click.Choice(["simple", "standard"]),
    default="simple",
    show_default=True,
    help="Validation depth (standard adds CheckM2 + GTDB-Tk).",
)
@click.option(
    "--output-dir",
    type=click.Path(file_okay=False),
    default=None,
    help="Directory for validation.json (default: <contigs-dir>/../taxonomy).",
)
def validate_cmd(contigs: str, mode: str, output_dir: str | None) -> None:
    """Validate an assembly with marker genes, optionally CheckM2/GTDB-Tk."""
    import json

    from .analysis.taxonomic_validator import validate_genome

    result = validate_genome(contigs, mode=mode, output_dir=output_dir)
    print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2))


# Register commands
main.add_command(species_cmd, name="species")
main.add_command(db_group, name="db")
main.add_command(validate_cmd, name="validate")

# Imported here, not at the top: gside.mcp pulls in mcp_schemas, which reads
# MODES/TIERS from this module, so they must already exist when it loads
# (unavoidable circular-import constraint).
from gside.mcp import register_mcp_command  # noqa: E402

register_mcp_command(main)


class _OutputFileWriter:
    """Chunk sink behind ``species --output`` (gapit OutputFileWriter semantics):
    opens truncating on the first chunk, and every chunk is flushed, so an
    in-batch failure still persists the already-streamed prefix."""

    def __init__(self, path: Path) -> None:
        self._path = path
        self._handle: TextIO | None = None

    def __call__(self, chunk: str) -> None:
        handle = self._handle or self._path.open("w", encoding="utf-8", newline="\n")
        self._handle = handle
        handle.write(chunk)
        handle.flush()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()


def _run_species(
    contigs_list: list[str],
    mode: str,
    db_dir: str | None,
    output_format: str = "json",
    output: str | None = None,
    jobs: int = 1,
    threads: int = 4,
    no_header: bool = False,
    auto_setup: bool = False,
    quiet: bool = False,
) -> int:
    from gside.progress import species_progress

    def identify(contigs: str) -> dict[str, Any]:
        payload = _identify_one(contigs, mode, db_dir, threads)
        if auto_setup:
            payload = _auto_setup_retry(contigs, payload, mode, db_dir, threads)
        return payload

    if jobs > 1:
        from concurrent.futures import ThreadPoolExecutor

        with (
            ThreadPoolExecutor(max_workers=jobs) as pool,
            species_progress(
                contigs_list, quiet=quiet, description=f"Identifying ({jobs} jobs)"
            ) as (bar, _tid),
        ):
            return _emit_species(
                _tracked(pool.map(identify, contigs_list), bar),
                len(contigs_list),
                output_format,
                output,
                no_header,
            )
    with species_progress(contigs_list, quiet=quiet) as (bar, _tid):
        return _emit_species(
            _tracked(map(identify, contigs_list), bar),
            len(contigs_list),
            output_format,
            output,
            no_header,
        )


def _tracked(source: Iterator[dict[str, Any]], bar: Any) -> Iterator[dict[str, Any]]:
    for payload in source:
        yield payload
        bar.advance()


def _emit_species(
    payload_source: Iterator[dict[str, Any]],
    n_inputs: int,
    output_format: str,
    output: str | None,
    no_header: bool = False,
) -> int:
    if output_format == "json":
        payloads = list(payload_source)
        report: Any = payloads[0] if n_inputs == 1 else payloads
        text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        if output is None:
            print(text, end="")
        else:
            writer = _OutputFileWriter(Path(output))
            writer(text)
            writer.close()
            print(f"Results written to {output}", file=sys.stderr)
        return 0

    markdown = output_format == "md"

    def chunks() -> Iterator[str]:
        if not (no_header and not markdown):
            yield _table_header(markdown) + "\n"
        for payload in payload_source:
            yield _render_row(payload, markdown) + "\n"

    if output is None:
        for chunk in chunks():
            print(chunk, end="", flush=True)
        return 0

    writer = _OutputFileWriter(Path(output))
    try:
        for chunk in chunks():
            writer(chunk)
    finally:
        writer.close()
    print(f"Results written to {output}", file=sys.stderr)
    return 0


TABLE_COLUMNS = [
    "contigs",
    "verdict_species",
    "verdict_confidence",
    "verdict_basis",
    "marker_species",
    "marker_confidence",
    "marker_rule",
    "panel_species",
    "panel_confidence",
    "mash_refseq_species",
    "mash_refseq_confidence",
    "sourmash_species",
    "sourmash_confidence",
    "errors",
]


def _method_cell(methods: dict[str, Any], mode: str, key: str) -> str:
    payload = methods.get(mode, {})
    res = payload.get("result", payload)
    value = res.get(key, "")
    return "" if value is None else str(value)


def _table_header(markdown: bool) -> str:
    if not markdown:
        return "\t".join(TABLE_COLUMNS)
    return "| " + " | ".join(TABLE_COLUMNS) + " |\n|" + "|".join(["---"] * len(TABLE_COLUMNS)) + "|"


def _render_row(payload: dict[str, Any], markdown: bool) -> str:
    methods = payload.get("methods", {})
    verdict = payload.get("verdict", {})
    errors = "; ".join(f"{m}: {p.get('error', '')}" for m, p in methods.items() if "error" in p)
    cells = [
        payload.get("contigs", ""),
        verdict.get("species", ""),
        verdict.get("confidence", ""),
        ",".join(verdict.get("basis", [])),
        _method_cell(methods, "marker", "species"),
        _method_cell(methods, "marker", "confidence"),
        methods.get("marker", {}).get("matched_rule", ""),
        _method_cell(methods, "panel", "species"),
        _method_cell(methods, "panel", "confidence"),
        _method_cell(methods, "mash_refseq", "species"),
        _method_cell(methods, "mash_refseq", "confidence"),
        _method_cell(methods, "sourmash", "species"),
        _method_cell(methods, "sourmash", "confidence"),
        errors,
    ]

    def _clean(cell: str) -> str:
        return cell.replace("|", "/").replace("\n", " ").replace("\t", " ").replace("\r", " ")

    clean = [_clean(c) for c in cells]
    return ("| " + " | ".join(clean) + " |") if markdown else "\t".join(clean)


def _render_table(payloads: list[dict[str, Any]], markdown: bool) -> str:
    lines = [_table_header(markdown)]
    lines += [_render_row(payload, markdown) for payload in payloads]
    return "\n".join(lines)


def _structured_missing_db_fields(error: str) -> dict[str, str]:
    """Agent-friendly fields for a missing-database error (empty dict otherwise)."""
    marker = "(run: gside db setup --tier "
    idx = error.find(marker)
    if idx == -1:
        return {}
    tier = error[idx + len(marker) :].split(")", 1)[0].strip()
    if not tier:
        return {}
    from gside.db import _TIERS  # lazy: db.py is only needed for this lookup

    fields = {
        "suggested_command": f"gside db setup --tier {tier}",
        "missing_tier": tier,
    }
    size = _TIERS.get(tier, {}).get("size")
    if size:
        fields["estimated_size"] = str(size)
    return fields


def _setup_missing_tier(tier: str) -> None:
    """Install one database tier with progress routed to stderr."""
    from gside.db import db_setup

    with contextlib.redirect_stdout(sys.stderr):
        db_setup(tier)


def _auto_setup_retry(
    contigs: str,
    payload: dict[str, Any],
    mode: str,
    db_dir: str | None,
    threads: int,
) -> dict[str, Any]:
    """--auto-setup: install each missing tier once, retry that method once."""
    methods = payload.get("methods", {})
    for m, method_payload in methods.items():
        tier = method_payload.get("missing_tier") if isinstance(method_payload, dict) else None
        if not tier:
            continue
        try:
            _setup_missing_tier(tier)
        except Exception as e:  # noqa: BLE001 — CLI boundary: setup failure keeps the error entry
            logger.warning("auto-setup of tier %s failed: %s", tier, e)
            continue
        retry = _identify_one(contigs, m, db_dir, threads)
        methods[m] = retry["methods"].get(m, method_payload)
    payload["verdict"] = _verdict_for(methods, mode)
    return payload


def _verdict_for(results: dict[str, Any], mode: str) -> dict[str, Any]:
    if mode == "all" and len(results) > 1:
        return _arbitrate(results)
    first = next(iter(results.values()), {})
    return {
        "species": first.get("species", first.get("result", {}).get("species", "Unknown")),
        "confidence": first.get("confidence", first.get("result", {}).get("confidence", "low")),
        "basis": list(results),
    }


def _identify_one(contigs: str, mode: str, db_dir: str | None, threads: int = 4) -> dict[str, Any]:
    results: dict[str, Any] = {}

    modes = ["marker", "panel", "mash_refseq", "sourmash", "gtdbtk"] if mode == "all" else [mode]
    for m in modes:
        try:
            if m == "marker":
                from gside.analysis.multigene_identifier import identify_multigene

                results[m] = identify_multigene(contigs, threads=threads).to_dict()
            elif m in ("panel", "mash_refseq"):
                from gside.analysis.ani_identifier import identify_by_ani

                results[m] = identify_by_ani(
                    contigs, mode=m, db_dir=db_dir, threads=threads
                ).to_dict()
            elif m == "sourmash":
                from gside.analysis.sourmash_identifier import (
                    identify_by_sourmash as identify_sourmash,
                )

                results[m] = identify_sourmash(contigs, db_dir=db_dir).to_dict()
            elif m == "gtdbtk":
                from gside.analysis.gtdbtk_identifier import identify_by_gtdbtk

                results[m] = identify_by_gtdbtk(contigs, threads=threads).to_dict()
        except Exception as e:  # noqa: BLE001 — CLI 边界：单法失败不终止其他方法
            entry = {"method": m, "error": str(e)[:300]}
            entry.update(_structured_missing_db_fields(entry["error"]))
            results[m] = entry

    verdict = _verdict_for(results, mode)

    payload = {
        "analysis_type": "species_identification",
        "tool": "gside",
        "version": __version__,
        "contigs": contigs,
        "methods": results,
        "verdict": verdict,
    }
    return payload


def _arbitrate(results: dict[str, Any]) -> dict[str, Any]:
    """Layer-first arbitration: gtdbtk (phylogenetic) > panel/mash/sourmash (ANI) > marker."""
    _LAYER = {"gtdbtk": 3, "panel": 2, "mash_refseq": 2, "sourmash": 2, "marker": 1}
    best = None
    for mode, payload in results.items():
        if "error" in payload:
            continue
        res = payload.get("result", payload)
        species = res.get("species", "Unknown")
        conf = res.get("confidence", "low")
        score = (_LAYER.get(mode, 0), 1 if conf == "high" else 0)
        if species != "Unknown" and (best is None or score > best[0]):
            best = (score, species, conf, mode)
    if best is None:
        return {"species": "Unknown", "confidence": "low", "basis": list(results)}
    return {"species": best[1], "confidence": best[2], "basis": [best[3]]}


if __name__ == "__main__":
    main()
