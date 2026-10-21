"""Command-line interface for gside."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import logging
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any, TextIO

import click

from gside import __version__
from gside.analysis.species_name import UNKNOWN, display_name, parse_species_name

HELP_SETTINGS = {"help_option_names": ["-h", "--help"]}

MODES = ["marker", "panel", "mash_refseq", "sourmash", "gtdbtk", "all"]
META_MODES = ["kraken2", "sylph", "all"]
TIERS = ["mini", "panel", "mash", "sourmash", "meta", "meta_sylph", "all"]


def _db_dir_option(func: Any) -> Any:
    """``-D/--db-dir`` — shared by every command that reads reference databases."""
    return click.option(
        "-D",
        "--db-dir",
        type=click.Path(file_okay=False),
        default=None,
        help="Database root (default: $GSIDE_DB_DIR or data/db).",
    )(func)


def _format_option(choices: list[str], help_text: str) -> Any:
    """``-f/--format`` factory so every command spells the flag identically."""

    def decorate(func: Any) -> Any:
        return click.option(
            "-f",
            "--format",
            "output_format",
            type=click.Choice(choices),
            default="json",
            show_default=True,
            help=help_text,
        )(func)

    return decorate


def _threads_option(help_text: str) -> Any:
    """``-t/--threads`` factory (per-tool wording differs, the bounds do not)."""

    def decorate(func: Any) -> Any:
        return click.option(
            "-t",
            "--threads",
            type=click.IntRange(min=1),
            default=4,
            show_default=True,
            help=help_text,
        )(func)

    return decorate


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
@_db_dir_option
@_format_option(["json", "tsv", "md"], "Output format (tsv/md render one row per input file).")
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
@_threads_option("Threads per aligner invocation (blastn -num_threads / skani -t / mash -p).")
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


@main.command("meta", context_settings=HELP_SETTINGS)
@click.argument("reads", nargs=-1, required=False, type=click.Path(dir_okay=False))
@click.option(
    "-m",
    "--mode",
    type=click.Choice(META_MODES),
    default="all",
    show_default=True,
    help="Profiling method (all runs kraken2 + sylph, then merges).",
)
@_db_dir_option
@_format_option(["json", "tsv", "csv"], "Output format (tsv/csv render one row per species).")
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help=(
        "Write the report to PATH instead of stdout (truncates; json is written"
        " once at the end; stdout then carries no data)"
    ),
)
@click.option(
    "-s",
    "--samples",
    type=click.Path(dir_okay=False),
    default=None,
    help=(
        "Batch manifest: sample_id<TAB>reads [<TAB>reads2] per line (header optional)."
        " Output is one JSON object per sample, in manifest order."
    ),
)
@click.option(
    "-j",
    "--jobs",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Number of manifest samples to profile concurrently (order is always preserved).",
)
@click.option(
    "-n",
    "--no-header",
    is_flag=True,
    default=False,
    help="Suppress the TSV/CSV header row (json keeps its structure).",
)
@click.option(
    "-a",
    "--auto-setup",
    is_flag=True,
    default=False,
    help="Install missing database tiers before profiling (progress on stderr).",
)
@_threads_option("Threads per tool invocation (kraken2 --threads / sylph sketch -t).")
def meta_cmd(
    reads: tuple[str, ...],
    mode: str,
    db_dir: str | None,
    output_format: str,
    output: str | None,
    samples: str | None,
    jobs: int,
    no_header: bool,
    auto_setup: bool,
    threads: int,
) -> None:
    """Metagenomic species profiling from FASTQ reads, contigs, or a samples manifest."""
    if samples and reads:
        raise click.UsageError("pass either READS or --samples, not both")
    if not samples and not reads:
        raise click.UsageError("missing input: pass READS or --samples")
    if samples:
        try:
            manifest = _read_samples(samples)
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        sys.exit(
            _run_meta_samples(
                manifest,
                mode,
                db_dir,
                output_format=output_format,
                output=output,
                threads=threads,
                jobs=jobs,
                no_header=no_header,
                auto_setup=auto_setup,
            )
        )
    if jobs != 1:
        raise click.UsageError("--jobs only applies to --samples batches")
    sys.exit(
        _run_meta(
            list(reads),
            mode,
            db_dir,
            output_format=output_format,
            output=output,
            threads=threads,
            no_header=no_header,
            auto_setup=auto_setup,
        )
    )


# Register commands
main.add_command(species_cmd, name="species")
main.add_command(db_group, name="db")
main.add_command(validate_cmd, name="validate")
main.add_command(meta_cmd, name="meta")

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


def _write_text_document(text: str, output: str | None) -> int:
    """Emit rendered rows to stdout or to the ``-o`` file (truncating)."""
    if output is None:
        print(text, end="")
        return 0
    writer = _OutputFileWriter(Path(output))
    try:
        writer(text)
    finally:
        writer.close()
    print(f"Results written to {output}", file=sys.stderr)
    return 0


def _write_json_document(payload: Any, output: str | None) -> int:
    """Emit one JSON document to stdout or to the ``-o`` file (truncating)."""
    return _write_text_document(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", output)


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
        return _write_json_document(report, output)

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
    info = _TIERS.get(tier, {})
    size = info.get("size")
    if size:
        fields["estimated_size"] = str(size)
    if info.get("help_url"):
        fields["help_url"] = str(info["help_url"])
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


def _auto_setup_retry_meta(
    files: list[str],
    payload: dict[str, Any],
    mode: str,
    db_dir: str | None,
    threads: int,
) -> dict[str, Any]:
    """--auto-setup for meta: install each missing tier once, retry that method once."""
    if db_dir is not None:
        # A custom database root is user-managed; never write outside it.
        return payload
    methods = payload.get("methods", {})
    input_type = payload.get("input", {}).get("type", _detect_meta_input(files))
    for m, entry in methods.items():
        tier = entry.get("missing_tier") if isinstance(entry, dict) else None
        if not tier:
            continue
        try:
            _setup_missing_tier(tier)
        except Exception as e:  # noqa: BLE001 — CLI boundary: setup failure keeps the error entry
            logger.warning("auto-setup of tier %s failed: %s", tier, e)
            continue
        methods[m] = _run_meta_method(m, files, db_dir, threads, input_type)
    payload["verdict"] = _meta_verdict(methods)
    return payload


def _species_of(payload: dict[str, Any]) -> tuple[str, str, dict[str, Any]]:
    """(species, confidence, result-dict) from a method payload."""
    res = payload.get("result", payload)
    return res.get("species", "Unknown"), res.get("confidence", "low"), res


def _verdict_payload(
    species: str,
    confidence: str,
    strain: Any,
    basis: list[str],
    supporting: list[str],
    agreement: str,
) -> dict[str, Any]:
    """Verdict object shape shared by single-mode and ``--mode all`` runs."""
    return {
        "species": species,
        "display_name": display_name(species),
        "strain": strain or None,
        "confidence": confidence,
        "basis": basis,
        "supporting": supporting,
        "agreement": agreement,
    }


def _verdict_for(results: dict[str, Any], mode: str) -> dict[str, Any]:
    if mode == "all" and len(results) > 1:
        return _arbitrate(results)
    first = next(iter(results.values()), {})
    if not first or "error" in first:
        return _verdict_payload("Unknown", "low", None, list(results), [], "0/1")
    species, conf, res = _species_of(first)
    return _verdict_payload(species, conf, res.get("strain"), list(results), [], "1/1")


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


_LAYER = {"gtdbtk": 3, "panel": 2, "mash_refseq": 2, "sourmash": 2, "marker": 1}


def _arbitrate(results: dict[str, Any]) -> dict[str, Any]:
    """Layer-first arbitration: gtdbtk (phylogenetic) > panel/mash/sourmash (ANI) > marker.

    Within a layer ``medium`` now outranks ``low``. ``basis`` keeps the single
    winning method (bacmap maps it to a validated ``method`` field), while
    ``supporting``/``agreement`` expose corroboration from the other tiers.
    """
    best_score: tuple[int, int] | None = None
    winner: tuple[str, str, str, dict[str, Any]] | None = None
    successful: dict[str, str] = {}

    for mode, payload in results.items():
        if "error" in payload:
            continue
        species, conf, res = _species_of(payload)
        successful[mode] = species
        if species == "Unknown":
            continue
        score = (_LAYER.get(mode, 0), _CONF_RANK.get(conf, 0))
        if best_score is None or score > best_score:
            best_score = score
            winner = (mode, species, conf, res)

    if winner is None:
        return _verdict_payload("Unknown", "low", None, list(results), [], f"0/{len(successful)}")

    mode, species, conf, res = winner
    supporting = sorted(m for m, s in successful.items() if m != mode and s == species)
    agreement = f"{1 + len(supporting)}/{len(successful)}"
    return _verdict_payload(species, conf, res.get("strain"), [mode], supporting, agreement)


_META_FASTA_SUFFIXES = (".fna", ".fa", ".fasta")
_META_TSV_HEADER = "method\tname\tabundance\tconfidence\tani"
_CONF_RANK = {"high": 2, "medium": 1, "low": 0}
_METHOD_RANK = {"sylph": 1, "kraken2": 0}


def _detect_meta_input(files: list[str]) -> str:
    """2 files = paired FASTQ; 1 FASTA-suffixed file = contigs; else single FASTQ."""
    if len(files) == 2:
        return "fastq_paired"
    if len(files) == 1 and Path(files[0]).suffix in _META_FASTA_SUFFIXES:
        return "fasta"
    return "fastq_single"


def _run_meta_method(
    mode: str,
    files: list[str],
    db_dir: str | None,
    threads: int,
    input_type: str,
) -> dict[str, Any]:
    from gside.analysis.sylph_identifier import FASTA_INPUT_ERROR

    if mode == "sylph" and input_type == "fasta":
        return {"method": "sylph", "error": FASTA_INPUT_ERROR}
    try:
        if mode == "kraken2":
            from gside.analysis.kraken2_identifier import identify_by_kraken2

            return identify_by_kraken2(files, db_dir=db_dir, threads=threads).to_dict()
        from gside.analysis.sylph_identifier import identify_by_sylph

        return identify_by_sylph(files, db_dir=db_dir, threads=threads).to_dict()
    except Exception as e:  # noqa: BLE001 — CLI boundary: method failure recorded, others continue
        entry = {"method": mode, "error": str(e)[:300]}
        entry.update(_structured_missing_db_fields(entry["error"]))
        return entry


def _meta_verdict(methods: dict[str, Any]) -> dict[str, Any]:
    """Merge species lists from successful methods into canonical rows.

    Rows are keyed by the canonical ``Genus_species`` label (the same contract
    as ``gside species``), so sylph's strain-level names and kraken2's
    taxonomy names collapse to one entry; on duplicates the higher-confidence
    method wins, with sylph (ANI) outranking kraken2 (fraction) at equal
    confidence. Species-group placeholders (``Campylobacter sp.``) are dropped.
    """
    merged: dict[str, dict[str, Any]] = {}
    basis: list[str] = []
    for mode, payload in methods.items():
        if "error" in payload:
            continue
        basis.append(mode)
        res = payload.get("result", payload)
        entries = res.get("top_species") or res.get("species") or []
        for entry in entries:
            raw_name = entry.get("name", "")
            if not raw_name:
                continue
            parsed = parse_species_name(raw_name)
            if parsed.label == UNKNOWN:
                continue
            conf = entry.get("confidence", "low")
            row = {
                "name": parsed.label,
                "display_name": parsed.display,
                "strain": parsed.strain or None,
                "method": mode,
                "abundance": entry.get("fraction", entry.get("abundance", 0.0)),
                "confidence": conf,
                "ani": entry.get("ani"),
            }
            if raw_name != parsed.label:
                row["species_raw"] = raw_name

            prev = merged.get(parsed.label)
            wins = prev is None or _CONF_RANK[conf] > _CONF_RANK[prev["confidence"]]
            if not wins and prev is not None and _CONF_RANK[conf] == _CONF_RANK[prev["confidence"]]:
                wins = _METHOD_RANK.get(mode, 0) > _METHOD_RANK.get(prev["method"], 0)
            if wins:
                merged[parsed.label] = row
    top = sorted(merged.values(), key=lambda r: r["abundance"], reverse=True)
    return {"top_species": top, "richness": len(top), "basis": basis}


_META_COLUMNS = ("method", "name", "abundance", "confidence", "ani")


def _meta_row_cells(payload: dict[str, Any]) -> list[list[str]]:
    """One row per merged species, values already stringified for tsv/csv."""
    rows: list[list[str]] = []
    for row in payload.get("verdict", {}).get("top_species", []):
        ani = "" if row.get("ani") is None else str(row["ani"])
        abundance = f"{row['abundance']:.6g}"
        rows.append([row["method"], row["name"], abundance, row["confidence"], ani])
    return rows


def _render_meta_table(
    payloads: list[dict[str, Any]], fmt: str, no_header: bool, with_sample: bool
) -> str:
    header = (["sample"] if with_sample else []) + list(_META_COLUMNS)
    if fmt == "csv":
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        if not no_header:
            writer.writerow(header)
        for payload in payloads:
            prefix = [payload.get("sample", "")] if with_sample else []
            for cells in _meta_row_cells(payload):
                writer.writerow(prefix + cells)
        return out.getvalue()

    lines = [] if no_header else ["\t".join(header)]
    for payload in payloads:
        prefix = [payload.get("sample", "")] if with_sample else []
        lines += ["\t".join(prefix + cells) for cells in _meta_row_cells(payload)]
    return "\n".join(lines) + "\n"


def _emit_meta(
    payload: dict[str, Any],
    output_format: str,
    output: str | None,
    no_header: bool = False,
) -> int:
    if output_format == "json":
        return _write_json_document(payload, output)
    text = _render_meta_table([payload], output_format, no_header, with_sample=False)
    return _write_text_document(text, output)


def _emit_meta_batch(
    payloads: list[dict[str, Any]],
    output_format: str,
    output: str | None,
    no_header: bool = False,
) -> int:
    if output_format == "json":
        return _write_json_document(payloads, output)
    text = _render_meta_table(payloads, output_format, no_header, with_sample=True)
    return _write_text_document(text, output)


def _meta_payload(
    files: list[str],
    mode: str,
    db_dir: str | None,
    threads: int = 4,
) -> dict[str, Any]:
    """Profile one sample into the JSON payload (shared by CLI, MCP, batch)."""
    input_type = _detect_meta_input(files)
    modes = ["kraken2", "sylph"] if mode == "all" else [mode]
    methods = {m: _run_meta_method(m, files, db_dir, threads, input_type) for m in modes}
    return {
        "analysis_type": "metagenomic_profiling",
        "tool": "gside",
        "version": __version__,
        "input": {"files": files, "type": input_type},
        "methods": methods,
        "verdict": _meta_verdict(methods),
    }


def _run_meta(
    files: list[str],
    mode: str,
    db_dir: str | None,
    output_format: str = "json",
    output: str | None = None,
    threads: int = 4,
    no_header: bool = False,
    auto_setup: bool = False,
) -> int:
    payload = _meta_payload(files, mode, db_dir, threads)
    if auto_setup:
        payload = _auto_setup_retry_meta(files, payload, mode, db_dir, threads)
    return _emit_meta(payload, output_format, output, no_header)


def _read_samples(path: str) -> list[dict[str, Any]]:
    """Parse a ``sample_id<TAB>r1[<TAB>r2]`` manifest (header optional)."""
    manifest = Path(path)
    if not manifest.is_file():
        raise ValueError(f"samples manifest not found: {manifest}")
    samples: list[dict[str, Any]] = []
    seen: set[str] = set()
    for lineno, line in enumerate(manifest.read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        cols = [c.strip() for c in line.split("\t")]
        if lineno == 1 and cols[0].lower() in ("sample", "sample_id"):
            continue
        if len(cols) < 2 or not cols[0] or not cols[1]:
            raise ValueError(f"{manifest}:{lineno}: expected sample_id<TAB>reads [<TAB>reads2]")
        if len(cols) > 3:
            raise ValueError(f"{manifest}:{lineno}: too many columns (max 3)")
        sample_id = cols[0]
        if sample_id in seen:
            raise ValueError(f"{manifest}:{lineno}: duplicate sample id {sample_id!r}")
        seen.add(sample_id)
        samples.append({"id": sample_id, "files": [c for c in cols[1:] if c]})
    if not samples:
        raise ValueError(f"{manifest}: no samples found")
    return samples


def _run_meta_samples(
    samples: list[dict[str, Any]],
    mode: str,
    db_dir: str | None,
    output_format: str = "json",
    output: str | None = None,
    threads: int = 4,
    jobs: int = 1,
    no_header: bool = False,
    auto_setup: bool = False,
) -> int:
    """Profile every manifest sample, preserving manifest order in the output."""

    def identify(sample: dict[str, Any]) -> dict[str, Any]:
        payload = _meta_payload(sample["files"], mode, db_dir, threads)
        if auto_setup:
            payload = _auto_setup_retry_meta(sample["files"], payload, mode, db_dir, threads)
        return {"sample": sample["id"], **payload}

    if auto_setup and samples:
        # Install tiers once, serially, before the parallel fan-out.
        first = identify(samples[0])
        rest = samples[1:]
        if jobs > 1 and rest:
            from concurrent.futures import ThreadPoolExecutor

            with ThreadPoolExecutor(max_workers=jobs) as pool:
                payloads = [first, *pool.map(identify, rest)]
        else:
            payloads = [first, *(identify(sample) for sample in rest)]
        return _emit_meta_batch(payloads, output_format, output, no_header)

    if jobs > 1:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=jobs) as pool:
            payloads = list(pool.map(identify, samples))
    else:
        payloads = [identify(sample) for sample in samples]
    return _emit_meta_batch(payloads, output_format, output, no_header)


if __name__ == "__main__":
    main()
