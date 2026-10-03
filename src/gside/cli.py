"""Command-line interface for gside."""

from __future__ import annotations

import json
import logging
import sys
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
)
@click.version_option(__version__, "--version", "-V")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
@click.option("--quiet", "-q", is_flag=True, help="Suppress non-error logging.")
@click.pass_context
def main(ctx: click.Context, verbose: bool, quiet: bool) -> None:
    """gside — Genome Species IDentification Engine.

    Multi-layer species identification from contigs. Prints a single
    JSON verdict per run to stdout.
    """
    if verbose and quiet:
        raise click.UsageError("--verbose and --quiet cannot be used together")
    setup_logging(verbose=verbose, quiet=quiet)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@main.command("species", context_settings=HELP_SETTINGS)
@click.argument("contigs", nargs=-1, required=True, type=click.Path(dir_okay=False))
@click.option(
    "--mode",
    type=click.Choice(MODES),
    default="marker",
    show_default=True,
    help="Identification method (all runs every method, then arbitrates).",
)
@click.option(
    "--db-dir",
    type=click.Path(file_okay=False),
    default=None,
    help="Database root (default: $GSIDE_DB_DIR or data/db).",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["json", "tsv", "md"]),
    default="json",
    show_default=True,
    help="Output format (tsv/md render one row per input file).",
)
def species_cmd(
    contigs: tuple[str, ...], mode: str, db_dir: str | None, output_format: str
) -> None:
    """Identify species from assembled contigs FASTA (accepts multiple files)."""
    sys.exit(_run_species(list(contigs), mode, db_dir, output_format))


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
def db_status_cmd() -> None:
    """Show database readiness status."""
    from .db import run_db_command

    sys.exit(run_db_command(["status"]))


@db_group.command("setup", context_settings=HELP_SETTINGS)
@click.option(
    "--tier",
    type=click.Choice(TIERS),
    default="panel",
    show_default=True,
    help="Database tier to install (all = panel + mash; sourmash is separate).",
)
@click.option(
    "--source",
    default="",
    help="Copy from an existing local database directory.",
)
def db_setup_cmd(tier: str, source: str) -> None:
    """Install or update databases."""
    from .db import run_db_command

    args = ["setup", tier]
    if source:
        args += ["--source", source]
    sys.exit(run_db_command(args))


@db_group.command("list", context_settings=HELP_SETTINGS)
def db_list_cmd() -> None:
    """List available database tiers."""
    from .db import run_db_command

    sys.exit(run_db_command(["list"]))


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
    contigs_list: list[str], mode: str, db_dir: str | None, output_format: str = "json"
) -> int:
    payloads = [_identify_one(contigs, mode, db_dir) for contigs in contigs_list]
    if output_format == "json":
        if len(payloads) == 1:
            print(json.dumps(payloads[0], ensure_ascii=False, indent=2))
        else:
            print(json.dumps(payloads, ensure_ascii=False, indent=2))
    else:
        print(_render_table(payloads, markdown=(output_format == "md")))
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
