"""Command-line interface for gside."""

from __future__ import annotations

import contextlib
import json
import logging
import sys
from pathlib import Path
from typing import Any

import click

from gside import __version__

HELP_SETTINGS = {"help_option_names": ["-h", "--help"]}

MODES = ["marker", "panel", "mash_refseq", "sourmash", "all"]
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
    help="Write the report to PATH instead of stdout (truncates; stdout then carries no data).",
)
def species_cmd(
    contigs: tuple[str, ...],
    mode: str,
    db_dir: str | None,
    output_format: str,
    output: str | None,
) -> None:
    """Identify species from assembled contigs FASTA (accepts multiple files)."""
    sys.exit(_run_species(list(contigs), mode, db_dir, output_format, output))


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
    type=click.Choice(["table", "tsv", "json"]),
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
    type=click.Choice(["table", "tsv", "json"]),
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


def _run_species(
    contigs_list: list[str],
    mode: str,
    db_dir: str | None,
    output_format: str = "json",
    output: str | None = None,
) -> int:
    payloads = [_identify_one(contigs, mode, db_dir) for contigs in contigs_list]
    if output_format == "json":
        report: Any = payloads[0] if len(payloads) == 1 else payloads
        text = json.dumps(report, ensure_ascii=False, indent=2)
    else:
        text = _render_table(payloads, markdown=(output_format == "md"))
    if output is None:
        print(text)
    else:
        Path(output).write_text(text + "\n", encoding="utf-8")
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


def _render_table(payloads: list[dict[str, Any]], markdown: bool) -> str:
    rows: list[list[str]] = []
    for payload in payloads:
        methods = payload.get("methods", {})
        verdict = payload.get("verdict", {})
        errors = "; ".join(f"{m}: {p.get('error', '')}" for m, p in methods.items() if "error" in p)
        rows.append(
            [
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
        )

    def _clean(cell: str) -> str:
        return cell.replace("|", "/").replace("\n", " ").replace("\t", " ").replace("\r", " ")

    clean = [[_clean(c) for c in row] for row in rows]
    if not markdown:
        return "\n".join(["\t".join(TABLE_COLUMNS)] + ["\t".join(row) for row in clean])
    lines = [
        "| " + " | ".join(TABLE_COLUMNS) + " |",
        "|" + "|".join(["---"] * len(TABLE_COLUMNS)) + "|",
    ]
    lines += ["| " + " | ".join(row) + " |" for row in clean]
    return "\n".join(lines)


def _identify_one(contigs: str, mode: str, db_dir: str | None) -> dict[str, Any]:
    results: dict[str, Any] = {}

    modes = ["marker", "panel", "mash_refseq", "sourmash"] if mode == "all" else [mode]
    for m in modes:
        try:
            if m == "marker":
                from gside.analysis.multigene_identifier import identify_multigene

                results[m] = identify_multigene(contigs).to_dict()
            elif m in ("panel", "mash_refseq"):
                from gside.analysis.ani_identifier import identify_by_ani

                results[m] = identify_by_ani(contigs, mode=m, db_dir=db_dir).to_dict()
            elif m == "sourmash":
                from gside.analysis.sourmash_identifier import (
                    identify_by_sourmash as identify_sourmash,
                )

                results[m] = identify_sourmash(contigs, db_dir=db_dir).to_dict()
        except Exception as e:  # noqa: BLE001 — CLI 边界：单法失败不终止其他方法
            results[m] = {"method": m, "error": str(e)[:300]}

    if mode == "all" and len(results) > 1:
        verdict = _arbitrate(results)
    else:
        first = next(iter(results.values()), {})
        verdict = {
            "species": first.get("species", first.get("result", {}).get("species", "Unknown")),
            "confidence": first.get("confidence", first.get("result", {}).get("confidence", "low")),
            "basis": list(results),
        }

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
    """ANI-layer (panel/mash/sourmash) high-confidence beats marker; else marker."""
    _LAYER = {"panel": 2, "mash_refseq": 2, "sourmash": 2, "marker": 1}
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
