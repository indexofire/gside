"""Command-line interface for gside."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import logging
import re
import sys
import threading
from collections.abc import Iterator
from enum import StrEnum
from pathlib import Path
from typing import Any, TextIO

import click

from gside import __version__
from gside.analysis.species_name import display_name

HELP_SETTINGS = {"help_option_names": ["-h", "--help"]}

_ISOLATE_MODES = ["marker", "ani", "mash", "sourmash", "gtdbtk"]
# Legacy mode spellings accepted at the CLI surface and normalized
# immediately; the canonical names mirror the tier directories.
MODE_ALIASES = {"panel": "ani", "mash_refseq": "mash"}
MODES = [*_ISOLATE_MODES, "all"]
SCAN_MODES = MODES
TIERS = ["mini", "panel", "mash", "sourmash", "all"]


def _db_dir_option(func: Any) -> Any:
    """``-D/--db-dir`` — shared by every command that reads reference databases."""
    return click.option(
        "-D",
        "--db-dir",
        type=click.Path(file_okay=False),
        default=None,
        help="Database root (default: $GSIDE_DB_DIR or data/db).",
    )(func)


def _marker_db_option(func: Any) -> Any:
    """``--marker-db`` — use a custom marker database directory."""
    return click.option(
        "--marker-db",
        type=click.Path(exists=True, file_okay=False),
        default=None,
        help="Custom marker database directory (markers.fasta + marker_rules.yaml).",
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


# --------------------------------------------------------------- scan modes

_FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".faa", ".ffn", ".frn")
_FASTQ_SUFFIXES = (".fq", ".fastq")


class InputKind(StrEnum):
    """Sequence input format inferred from a file extension."""

    FASTA = "fasta"
    FASTQ = "fastq"
    UNKNOWN = "unknown"


# Which input kinds each mode can consume. `all` means every species tier,
# so it needs FASTA.
_MODE_INPUTS: dict[str, frozenset[InputKind]] = {
    "marker": frozenset({InputKind.FASTA, InputKind.FASTQ}),
    "ani": frozenset({InputKind.FASTA}),
    "mash": frozenset({InputKind.FASTA, InputKind.FASTQ}),
    "sourmash": frozenset({InputKind.FASTA, InputKind.FASTQ}),
    "gtdbtk": frozenset({InputKind.FASTA}),
}


def _input_kind(path: str) -> InputKind:
    """'fasta' / 'fastq' / 'unknown' from the extension (``.gz`` aware)."""
    name = path.lower()
    for suffix in (".gz", ".bz2"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if name.endswith(_FASTA_SUFFIXES):
        return InputKind.FASTA
    if name.endswith(_FASTQ_SUFFIXES):
        return InputKind.FASTQ
    return InputKind.UNKNOWN


_MATE_RE = re.compile(r"^(.*?)[._-](?:R)?([12])(?:_001)?$", re.IGNORECASE)


def _mate_key(path: str) -> tuple[str, int] | None:
    """`S_R1.fq.gz` → ("S", 1); None when the name has no mate marker."""
    name = Path(path).name
    for suffix in (".gz", ".bz2"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in (".fastq", ".fq"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    match = _MATE_RE.match(name)
    if not match:
        return None
    return match.group(1), int(match.group(2))


def _group_read_samples(files: list[str]) -> list[list[str]]:
    """Group FASTQ mates (`_R1/_R2`, `_1/_2`, `.1/.2`) into one sample.

    Order follows the first occurrence of each sample; mate 1 precedes mate 2
    in the produced pair. Files without a mate marker, FASTA files and
    unmatched mates each become their own single-file sample.
    """
    keyed = [(_mate_key(f) if _input_kind(f) == "fastq" else None, f) for f in files]
    bases: dict[str, dict[int, str]] = {}
    for key, path in keyed:
        if key:
            bases.setdefault(key[0].lower(), {})[key[1]] = path

    samples: list[list[str]] = []
    emitted: set[str] = set()
    for key, path in keyed:
        if key is None:
            samples.append([path])
            continue
        base = key[0].lower()
        if base in emitted:
            continue
        mates = bases[base]
        if 1 in mates and 2 in mates:
            samples.append([mates[1], mates[2]])
        else:
            samples.append([mates[1] if 1 in mates else mates[2]])
        emitted.add(base)
    return samples


def _supported_inputs(mode: str) -> frozenset[InputKind]:
    return _MODE_INPUTS[mode]


def _split_modes(raw: str) -> list[str]:
    """``marker+sourmash`` → ['marker', 'sourmash'] (order kept, dupes dropped).

    Legacy spellings (``panel``, ``mash_refseq``) normalize to the canonical
    ``ani``/``mash`` before any validation, so the rest of the pipeline only
    ever sees canonical mode names.
    """
    parts = (
        MODE_ALIASES.get(part.strip(), part.strip()) for part in raw.split("+") if part.strip()
    )
    return list(dict.fromkeys(parts))


def _expand_modes(raw: str) -> list[str]:
    parts = _split_modes(raw)
    if parts == ["all"]:
        return list(_ISOLATE_MODES)
    return parts


def _resolve_modes(raw: str) -> list[str]:
    """Validated expansion used before any tool runs."""
    parts = _split_modes(raw)
    if not parts:
        raise click.UsageError("--mode is empty")
    if "all" in parts and len(parts) > 1:
        raise click.UsageError("'all' cannot be combined with other modes")
    return _expand_modes(raw)


def _validate_scan(files: list[str], raw_mode: str) -> None:
    """Reject mode/input mismatches before any tool runs."""
    modes = _resolve_modes(raw_mode)
    for mode in modes:
        if mode not in _ISOLATE_MODES and mode != "all":
            raise click.UsageError(
                f"unknown mode {mode!r}; allowed: {', '.join([*_ISOLATE_MODES, 'all'])}"
            )

    kinds = {_input_kind(f) for f in files}
    if "unknown" in kinds:
        bad = next(f for f in files if _input_kind(f) == "unknown")
        raise click.UsageError(
            f"cannot infer input type from {bad!r}; "
            "expected FASTA (.fna/.fa/.fasta) or FASTQ (.fq/.fastq)"
        )

    for mode in modes:
        allowed = _supported_inputs(mode)
        unsupported = kinds - allowed
        if not unsupported:
            continue
        kind = sorted(unsupported)[0]
        if raw_mode == "all":
            hint = (
                "; all includes panel/gtdbtk, which need FASTA;"
                " for reads use e.g. -m marker+sourmash"
            )
        else:
            hint = ""
        raise click.UsageError(
            f"mode {mode!r} does not accept {kind.upper()} input"
            f" (supported: {', '.join(sorted(k.upper() for k in allowed))}){hint}"
        )


def _scan_dispatch(
    files: tuple[str, ...],
    mode: str | None,
    db_dir: str | None,
    marker_db: str | None,
    db_selection: tuple[str, ...],
    output_format: str,
    output: str | None,
    jobs: int,
    threads: int,
    no_header: bool,
    auto_setup: bool,
) -> None:
    paths = list(files)
    resolved = mode or "marker"
    # Normalize legacy spellings (panel→ani, mash_refseq→mash) once, here.
    resolved = "+".join(_split_modes(resolved))

    if not paths:
        raise click.UsageError("missing input: pass one or more FASTA/FASTQ files")
    grouped = _group_read_samples(paths)
    for sample in grouped:
        _validate_scan(sample, resolved)
    db_selections = _resolve_db_selections(db_selection, marker_db, resolved, db_dir)
    sys.exit(
        _run_species(
            grouped,
            resolved,
            db_dir,
            output_format,
            output,
            jobs,
            threads,
            no_header,
            auto_setup,
            marker_db=marker_db,
            db_selections=db_selections,
        )
    )


def _resolve_db_selections(
    db_selection: tuple[str, ...],
    marker_db: str | None,
    resolved_mode: str,
    db_dir: str | None,
) -> dict[str, list[Path]]:
    """Fold ``-d/--db`` values (and the ``--marker-db`` alias) per mode.

    Returns ``{mode: [database paths]}`` for the modes with explicit
    selections; modes without one keep their default tier. Resolution and
    validation errors surface as CLI usage errors.
    """
    from gside.db_registry import resolve_selections

    values = list(db_selection)
    if marker_db is not None:
        values.append(f"marker:{marker_db}")
    if not values:
        return {}
    modes = _expand_modes(resolved_mode)
    try:
        selections = resolve_selections(values, modes, root=db_dir)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    return {mode: [choice.path for choice in choices] for mode, choices in selections.items()}


def _scan_options(*, modes: list[str]) -> Any:
    """Shared option set for ``scan`` and its compatibility alias."""

    def decorate(func: Any) -> Any:
        decorators: list[Any] = [
            click.argument("files", nargs=-1, required=False, type=click.Path(dir_okay=False)),
            click.option(
                "-m",
                "--mode",
                default=None,
                metavar="MODE[+MODE...]",
                help=(f"Method or '+'-joined combination ({', '.join(modes)}); default: marker."),
            ),
            _db_dir_option,
            _marker_db_option,
            click.option(
                "-d",
                "--db",
                "db_selection",
                multiple=True,
                metavar="[MODE:]NAME",
                help=(
                    "Database(s) per mode: a name under the mode's tier directory"
                    " (data/db/D1_marker|D2_ani|D3_mash|D4_sourmash/<name>),"
                    " 'default', an existing directory path, or MODE:NAME when"
                    " several database-capable modes run together (prefix takes"
                    " the mode or its alias: panel:==ani:, mash_refseq:==mash:)."
                    " Repeatable; selections run and their results merge."
                ),
            ),
            _format_option(["json", "tsv", "md", "csv"], "Output format (rows for tsv/md/csv)."),
            click.option(
                "-o",
                "--output",
                type=click.Path(dir_okay=False),
                default=None,
                help=(
                    "Write the report to PATH instead of stdout (truncates; row formats"
                    " stream as inputs complete; json is written once at the end)"
                ),
            ),
            click.option(
                "-j",
                "--jobs",
                type=click.IntRange(min=1),
                default=1,
                show_default=True,
                help="Concurrent inputs (isolate batches).",
            ),
            _threads_option("Threads per tool invocation."),
            click.option(
                "-n",
                "--no-header",
                is_flag=True,
                default=False,
                help="Suppress the header row (tsv/csv/md).",
            ),
            click.option(
                "-a",
                "--auto-setup",
                is_flag=True,
                default=False,
                help="Install a missing database tier, then retry that method once.",
            ),
        ]
        for decorator in reversed(decorators):
            func = decorator(func)
        return func

    return decorate


@main.command("scan", context_settings=HELP_SETTINGS)
@_scan_options(modes=SCAN_MODES)
def scan_cmd(
    files: tuple[str, ...],
    mode: str | None,
    db_dir: str | None,
    marker_db: str | None,
    db_selection: tuple[str, ...],
    output_format: str,
    output: str | None,
    jobs: int,
    threads: int,
    no_header: bool,
    auto_setup: bool,
) -> None:
    """Identify the species of one or more isolate genomes.

    The input type is inferred from the extension and must match the mode
    (FASTA: marker/ani/mash/sourmash/gtdbtk/all; FASTQ:
    marker/mash/sourmash). Every mode can run several databases
    (``-d``; the tier directories D1_marker/D2_ani/D3_mash/D4_sourmash hold
    one subdirectory per database, ``default`` plus named ones). Metagenomic
    community profiling moved to the gmeta tool.
    """
    _scan_dispatch(
        files,
        mode,
        db_dir,
        marker_db,
        db_selection,
        output_format,
        output,
        jobs,
        threads,
        no_header,
        auto_setup,
    )


@main.command("species", hidden=True, context_settings=HELP_SETTINGS)
@_scan_options(modes=MODES)
def species_alias(
    files: tuple[str, ...],
    mode: str | None,
    db_dir: str | None,
    marker_db: str | None,
    db_selection: tuple[str, ...],
    output_format: str,
    output: str | None,
    jobs: int,
    threads: int,
    no_header: bool,
    auto_setup: bool,
) -> None:
    """Compatibility alias for ``gside scan`` (single organism)."""
    _scan_dispatch(
        files,
        mode,
        db_dir,
        marker_db,
        db_selection,
        output_format,
        output,
        jobs,
        threads,
        no_header,
        auto_setup,
    )


@click.group(
    "db",
    context_settings=HELP_SETTINGS,
    invoke_without_command=True,
)
@click.pass_context
def db_group(ctx: click.Context) -> None:
    """Manage reference databases (status/setup/list/validate)."""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        ctx.exit(0)


@db_group.command("status", hidden=True, context_settings=HELP_SETTINGS)
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
    """Hidden alias of ``db list`` (same mode → database view)."""
    from .db import format_status

    print(format_status(output_format))


@db_group.command("setup", context_settings=HELP_SETTINGS)
@click.option(
    "-t",
    "--tier",
    type=click.Choice(TIERS),
    default=None,
    help="Install by tier name (mini/panel/mash/sourmash/all). Legacy flag; prefer -d.",
)
@click.option(
    "-d",
    "--db",
    "db_selection",
    default=None,
    metavar="[MODE:]default",
    help="Install the default database of a mode: marker, panel, mash_refseq,"
    " sourmash (or their aliases ani/mash), MODE:default, or 'all'.",
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
def db_setup_cmd(
    tier: str | None,
    db_selection: str | None,
    source: str,
    output_format: str,
) -> None:
    """Install or update the default database of a mode (its tier)."""
    from .db import _is_setup_error, db_setup, run_db_command

    resolved = _resolve_setup_target(tier, db_selection)
    if isinstance(resolved, click.UsageError):
        raise resolved
    if output_format == "json":
        with contextlib.redirect_stdout(sys.stderr):
            results = db_setup(resolved, source)
        print(json.dumps(results, indent=2, ensure_ascii=False))
        sys.exit(1 if any(_is_setup_error(str(msg)) for msg in results.values()) else 0)

    args = ["setup", "--tier", resolved]
    if source:
        args += ["--source", source]
    sys.exit(run_db_command(args))


_SETUP_MODE_TO_TIER = {
    "marker": "mini",
    "ani": "panel",
    "mash": "mash",
    "sourmash": "sourmash",
    "all": "all",
    # legacy spellings
    "panel": "panel",
    "mash_refseq": "mash",
}


def _resolve_setup_target(tier: str | None, db_selection: str | None) -> str | click.UsageError:
    """``-d`` mode/alias/MODE:default vocabulary → tier name; ``-t`` legacy."""
    if db_selection is not None and tier is not None:
        return click.UsageError("-t/--tier and -d/--db are mutually exclusive")
    if db_selection is not None:
        raw = db_selection.strip()
        if raw.endswith(":default"):
            raw = raw[: -len(":default")]
        if raw == "default":
            return click.UsageError(
                "bare 'default' is ambiguous; name the mode: -d marker, -d panel,"
                " -d mash_refseq, -d sourmash"
            )
        mapped = _SETUP_MODE_TO_TIER.get(raw)
        if mapped is None:
            return click.UsageError(
                f"unknown database {db_selection!r}; -d installs a mode's default:"
                " marker, ani, mash, sourmash, all"
            )
        return mapped
    return tier if tier is not None else "panel"


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
    """List databases per mode (default + named) with readiness."""
    from .db import format_databases

    print(format_databases(output_format))


@db_group.command("validate", context_settings=HELP_SETTINGS)
@click.argument(
    "path",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["table", "json"]),
    default="table",
    show_default=True,
    help="Output format.",
)
def db_validate_cmd(path: Path, output_format: str) -> None:
    """Validate a custom marker database."""
    from gside.marker_db import MarkerDatabase

    try:
        database = MarkerDatabase.from_path(path)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    errors = database.validate()
    payload = {**database.to_dict(), "errors": errors}
    if output_format == "json":
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        print(f"Marker database: {database.name}")
        print(f"  path:    {database.root}")
        print(f"  version: {database.version}")
        if database.description:
            print(f"  about:   {database.description}")
        if errors:
            print("  errors:")
            for error in errors:
                print(f"    - {error}")
        else:
            print("  status:  valid")
    sys.exit(1 if errors else 0)


@main.command("validate", context_settings=HELP_SETTINGS)
@click.argument("contigs", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory for validation.json (created if needed). Required.",
)
@click.option(
    "-t",
    "--threads",
    type=click.IntRange(min=1),
    default=4,
    show_default=True,
    help="Threads for CheckM2 and GTDB-Tk.",
)
def validate_cmd(contigs: str, output_dir: str, threads: int) -> None:
    """Assess assembly quality and taxonomy (CheckM2 + GTDB-Tk)."""
    from .analysis.taxonomic_validator import validate_genome

    result = validate_genome(contigs, output_dir, threads=threads)
    print(json.dumps(result.to_dict(), ensure_ascii=False, indent=2))
    sys.exit(0 if result.status == "ok" else 1)


# Register commands
main.add_command(scan_cmd, name="scan")
main.add_command(species_alias, name="species")
main.add_command(db_group, name="db")
main.add_command(validate_cmd, name="validate")

# Imported here, not at the top: gside.mcp pulls in mcp_schemas, which reads
# MODES/TIERS from this module, so they must already exist when it loads
# (unavoidable circular-import constraint).
from gside.mcp import register_mcp_command  # noqa: E402

register_mcp_command(main)


class _WriterState(StrEnum):
    """Explicit file-sink lifecycle: IDLE -> OPEN -> CLOSED."""

    IDLE = "idle"
    OPEN = "open"
    CLOSED = "closed"


class _OutputFileWriter:
    """Chunk sink behind ``species --output`` (gspec OutputFileWriter semantics):
    opens truncating on the first chunk, and every chunk is flushed, so an
    in-batch failure still persists the already-streamed prefix."""

    def __init__(self, path: Path) -> None:
        self._path = path
        self._handle: TextIO | None = None
        self._state = _WriterState.IDLE

    @property
    def state(self) -> _WriterState:
        return self._state

    def __enter__(self) -> _OutputFileWriter:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def __call__(self, chunk: str) -> None:
        if self._state is _WriterState.CLOSED:
            raise RuntimeError(f"cannot write to closed output file: {self._path}")
        if self._handle is None:
            self._handle = self._path.open("w", encoding="utf-8", newline="\n")
            self._state = _WriterState.OPEN
        self._handle.write(chunk)
        self._handle.flush()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        self._state = _WriterState.CLOSED


def _write_text_document(text: str, output: str | None, bar: Any = None) -> int:
    """Emit rendered rows to stdout or to the ``-o`` file (truncating)."""
    if output is None:
        if bar is not None:
            bar.write(text)
        else:
            print(text, end="")
        return 0
    with _OutputFileWriter(Path(output)) as writer:
        writer(text)
    notice = f"Results written to {output}"
    bar.notice(notice) if bar is not None else print(notice, file=sys.stderr)
    return 0


def _write_json_document(payload: Any, output: str | None, bar: Any = None) -> int:
    """Emit one JSON document to stdout or to the ``-o`` file (truncating)."""
    text = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    return _write_text_document(text, output, bar)


def _run_species(
    contigs_list: list[str] | list[list[str]],
    mode: str,
    db_dir: str | None,
    output_format: str = "json",
    output: str | None = None,
    jobs: int = 1,
    threads: int = 4,
    no_header: bool = False,
    auto_setup: bool = False,
    quiet: bool = False,
    marker_db: str | None = None,
    db_selections: dict[str, list[Path]] | None = None,
) -> int:
    from gside.progress import species_progress

    samples = [list(item) if isinstance(item, (list, tuple)) else [item] for item in contigs_list]
    coordinator = _AutoSetupCoordinator() if auto_setup else None

    def make_identifier(bar: Any) -> Any:
        def identify(sample: list[str]) -> dict[str, Any]:
            target: str | list[str] = sample[0] if len(sample) == 1 else sample
            payload = _identify_one(target, mode, db_dir, threads, marker_db, db_selections)
            if coordinator is not None:
                payload = _auto_setup_retry(
                    target,
                    payload,
                    mode,
                    db_dir,
                    threads,
                    marker_db,
                    coordinator=coordinator,
                    bar=bar,
                    db_selections=db_selections,
                )
            return payload

        return identify

    if jobs > 1:
        from concurrent.futures import ThreadPoolExecutor

        with (
            ThreadPoolExecutor(max_workers=jobs) as pool,
            species_progress(samples, quiet=quiet, description=f"Identifying ({jobs} jobs)") as (
                bar,
                _tid,
            ),
        ):
            return _emit_species(
                _tracked(pool.map(make_identifier(bar), samples), bar),
                len(samples),
                output_format,
                output,
                no_header,
                bar,
            )
    with species_progress(samples, quiet=quiet) as (bar, _tid):
        return _emit_species(
            _tracked(map(make_identifier(bar), samples), bar),
            len(samples),
            output_format,
            output,
            no_header,
            bar,
        )


def _tracked(source: Iterator[dict[str, Any]], bar: Any) -> Iterator[dict[str, Any]]:
    for payload in source:
        yield payload
        bar.advance()


def _emit_chunks(chunks: Iterator[str], output: str | None, bar: Any = None) -> int:
    """Stream rendered chunks to stdout (flushed per row) or to ``-o``.

    When stdout is the destination and a live progress bar exists, every
    chunk goes through ``bar.write`` so the bar is suspended first; otherwise
    Rich's next refresh would erase the row.
    """
    if output is None:
        write = bar.write if bar is not None else lambda text: print(text, end="", flush=True)
        for chunk in chunks:
            write(chunk)
        return 0
    with _OutputFileWriter(Path(output)) as writer:
        for chunk in chunks:
            writer(chunk)
    notice = f"Results written to {output}"
    bar.notice(notice) if bar is not None else print(notice, file=sys.stderr)
    return 0


def _csv_line(cells: list[str]) -> str:
    buf = io.StringIO()
    csv.writer(buf, lineterminator="\n").writerow(cells)
    return buf.getvalue()


def _json_stream(payload_source: Iterator[dict[str, Any]], single: bool) -> Iterator[str]:
    """Stream a JSON document incrementally, one payload at a time.

    A single input yields one object; a batch yields ``[obj, obj, ...]`` with
    each object emitted as soon as it completes, so ``-o`` can be streamed too.
    """
    if single:
        first = next(payload_source, None)
        if first is not None:
            yield json.dumps(first, ensure_ascii=False, indent=2) + "\n"
        return
    yield "[\n"
    for index, payload in enumerate(payload_source):
        text = json.dumps(payload, ensure_ascii=False, indent=2)
        yield ("" if index == 0 else ",\n") + text + "\n"
    yield "]\n"


def _emit_species(
    payload_source: Iterator[dict[str, Any]],
    n_inputs: int,
    output_format: str,
    output: str | None,
    no_header: bool = False,
    bar: Any = None,
) -> int:
    if output_format == "json":
        return _emit_chunks(_json_stream(payload_source, n_inputs == 1), output, bar)

    markdown = output_format == "md"
    is_csv = output_format == "csv"

    def chunks() -> Iterator[str]:
        if not (no_header and not markdown):
            if is_csv:
                yield _csv_line(TABLE_COLUMNS)
            else:
                yield _table_header(markdown) + "\n"
        for payload in payload_source:
            if is_csv:
                yield _csv_line(_render_row(payload, markdown=False).split("\t"))
            else:
                yield _render_row(payload, markdown) + "\n"

    return _emit_chunks(chunks(), output, bar)


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
        _method_cell(methods, "ani", "species"),
        _method_cell(methods, "ani", "confidence"),
        _method_cell(methods, "mash", "species"),
        _method_cell(methods, "mash", "confidence"),
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


class _TierProvisionState(StrEnum):
    """Per-tier auto-setup lifecycle within a batch run."""

    PENDING = "pending"
    INSTALLING = "installing"
    INSTALLED = "installed"
    FAILED = "failed"


class _AutoSetupCoordinator:
    """Thread-safe batch state machine for ``--auto-setup``.

    Guarantees each missing tier is installed at most once per batch even
    under ``-j > 1``: concurrent workers wait while the first worker installs,
    a failed tier is never retried across remaining samples, and any live
    progress bar is suspended while setup progress streams to stderr.
    """

    def __init__(self) -> None:
        self._cond = threading.Condition()
        self._states: dict[str, _TierProvisionState] = {}

    def state_of(self, tier: str) -> _TierProvisionState:
        with self._cond:
            return self._states.get(tier, _TierProvisionState.PENDING)

    def provision(self, tier: str, bar: Any = None) -> bool:
        """Ensure *tier* was installed; returns True when retry should run."""
        with self._cond:
            while True:
                state = self._states.get(tier, _TierProvisionState.PENDING)
                if state is _TierProvisionState.INSTALLED:
                    return True
                if state is _TierProvisionState.FAILED:
                    return False
                if state is _TierProvisionState.INSTALLING:
                    self._cond.wait()
                    continue
                self._states[tier] = _TierProvisionState.INSTALLING
                break

        ok = False
        try:
            suspend = bar.suspend() if bar is not None else contextlib.nullcontext()
            with suspend:
                setup_result = _setup_missing_tier(tier)
            if isinstance(setup_result, dict):
                from gside.db import _is_setup_error

                ok = not any(_is_setup_error(str(msg)) for msg in setup_result.values())
            else:
                ok = True
        except Exception as e:  # noqa: BLE001 — CLI boundary: setup failure keeps the error entry
            logger.warning("auto-setup of tier %s failed: %s", tier, e)
            ok = False

        with self._cond:
            self._states[tier] = _TierProvisionState.INSTALLED if ok else _TierProvisionState.FAILED
            self._cond.notify_all()
        return ok


def _setup_missing_tier(tier: str) -> dict[str, str] | None:
    """Install one database tier with progress routed to stderr."""
    from gside.db import db_setup

    with contextlib.redirect_stdout(sys.stderr):
        return db_setup(tier)


def _auto_setup_retry(
    contigs: str | list[str],
    payload: dict[str, Any],
    mode: str,
    db_dir: str | None,
    threads: int,
    marker_db: str | None = None,
    *,
    coordinator: _AutoSetupCoordinator | None = None,
    bar: Any = None,
    db_selections: dict[str, list[Path]] | None = None,
) -> dict[str, Any]:
    """--auto-setup: install each missing tier once, retry that method once."""
    coord = coordinator if coordinator is not None else _AutoSetupCoordinator()
    methods = payload.get("methods", {})
    for m, method_payload in methods.items():
        tier = method_payload.get("missing_tier") if isinstance(method_payload, dict) else None
        if not tier:
            continue
        if not coord.provision(tier, bar=bar):
            continue
        if db_selections:
            retry = _identify_one(
                contigs, m, db_dir, threads, marker_db, db_selections=db_selections
            )
        elif marker_db is not None:
            retry = _identify_one(contigs, m, db_dir, threads, marker_db)
        else:
            retry = _identify_one(contigs, m, db_dir, threads)
        methods[m] = retry["methods"].get(m, method_payload)
    payload["verdict"] = _verdict_for(methods)
    return payload


def _species_of(payload: dict[str, Any]) -> tuple[str, str, dict[str, Any]]:
    """(species, confidence, result-dict) from a method payload."""
    res = payload.get("result", payload)
    return res.get("species", "Unknown"), res.get("confidence", "low"), res


# Deliberate genus-level labels: the evidence below species rank belongs to the
# typing tool (gspec), so the verdict says what it actually knows.
# Salmonella is NOT here anymore: the two species have their own marker rules
# (Salmonella_enterica / Salmonella_bongori); the bare label remains only as
# the genus-fallback rule output.
_GENUS_SCOPE_LABELS = frozenset({"Salmonella", "Shigella"})

# Species-complex labels: chromosomally inseparable groups where the honest
# call stops at complex level (the cereus group - B. anthracis vs B. cereus
# vs B. thuringiensis share >99% chromosomal identity; only the anthracis
# clade is separable via the lambdaBa03 marker).
_COMPLEX_SCOPE_LABELS = frozenset({"Bacillus_cereus_group"})

# Benchmark-informed close-relative groups. When successful methods disagree
# inside one group, the winner stays in ``species`` for backward compatibility,
# but the verdict gains ``scope: complex``, the competing labels in
# ``alternatives``, and a machine-readable ``ambiguity_reason``.
_TAXONOMIC_COMPLEXES: tuple[tuple[frozenset[str], str], ...] = (
    (
        frozenset(
            {
                "Escherichia_coli",
                "Shigella",
                "Shigella_flexneri",
                "Shigella_sonnei",
                "Shigella_boydii",
                "Shigella_dysenteriae",
            }
        ),
        "same_genomic_species_complex",
    ),
    (
        frozenset(
            {
                "Clostridium_botulinum",
                "Clostridium_haemolyticum",
                "Clostridium_combesii",
                "Clostridium_caccae",
                "Clostridium_novyi",
                "Clostridium_sporogenes",
            }
        ),
        "clostridium_botulinum_complex",
    ),
    (frozenset({"Enterococcus_faecium", "Enterococcus_lactis"}), "close_species_complex"),
    (
        frozenset({"Burkholderia_pseudomallei", "Burkholderia_mallei"}),
        "close_species_complex",
    ),
    (
        frozenset({"Vibrio_cholerae", "Vibrio_paracholerae"}),
        "close_species_complex",
    ),
    (
        frozenset({"Aeromonas_hydrophila", "Aeromonas_oralensis"}),
        "close_species_complex",
    ),
    (
        frozenset({"Legionella_pneumophila", "Legionella_fraseri", "Legionella_pascullei"}),
        "close_species_complex",
    ),
    (
        frozenset(
            {"Haemophilus_influenzae", "Haemophilus_parainfluenzae", "Haemophilus_aegyptius"}
        ),
        "close_species_complex",
    ),
    (
        frozenset(
            {
                "Bacillus_anthracis",
                "Bacillus_cereus",
                "Bacillus_cereus_group",
                "Bacillus_thuringiensis",
                "Bacillus_mycoides",
                "Bacillus_weihenstephanensis",
            }
        ),
        "bacillus_cereus_group_complex",
    ),
)

_VERDICT_EVIDENCE_KEYS = (
    "ani",
    "identity",
    "aligned_fraction",
    "f_unique_weighted",
)


class VerdictScope(StrEnum):
    """Evidence granularity of an arbitrated species verdict."""

    SPECIES = "species"
    GENUS = "genus"
    COMPLEX = "complex"
    UNKNOWN = "unknown"


def _scope_for_species(species: str) -> VerdictScope:
    """Evidence granularity of a verdict label."""
    if species == "Unknown":
        return VerdictScope.UNKNOWN
    if species in _GENUS_SCOPE_LABELS:
        return VerdictScope.GENUS
    if species in _COMPLEX_SCOPE_LABELS:
        return VerdictScope.COMPLEX
    return VerdictScope.SPECIES


def _complex_reason(labels: set[str]) -> str | None:
    """Ambiguity reason when every label belongs to one close-relative group."""
    for members, reason in _TAXONOMIC_COMPLEXES:
        if labels <= members:
            return reason
    return None


def _method_evidence(result: dict[str, Any]) -> dict[str, Any]:
    evidence = {key: result[key] for key in _VERDICT_EVIDENCE_KEYS if result.get(key) is not None}
    return evidence


def _verdict_semantics(
    species: str,
    successful: dict[str, str],
    results: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """``scope``/``alternatives``/``ambiguity_reason`` for an arbitrated call.

    Exact agreement keeps ``scope`` at species/genus level. A disagreement inside
    a known close-relative group upgrades scope to ``complex`` and records the
    competing labels with their methods and numeric evidence.
    """
    labels = {species, *(value for value in successful.values() if value != "Unknown")}
    reason = _complex_reason(labels) if len(labels) > 1 else None
    if reason is None:
        return {
            "scope": _scope_for_species(species),
            "alternatives": [],
            "ambiguity_reason": None,
        }

    grouped = sorted(label for label in labels if label != species)
    alternatives: list[dict[str, Any]] = []
    for alternative in grouped:
        modes = sorted(mode for mode, label in successful.items() if label == alternative)
        best_confidence = "low"
        evidence: dict[str, Any] = {}
        for mode in modes:
            _, confidence, result = _species_of(results[mode])
            if _CONF_RANK.get(confidence, 0) > _CONF_RANK.get(best_confidence, 0):
                best_confidence = confidence
            evidence[mode] = _method_evidence(result)
        alternatives.append(
            {
                "species": alternative,
                "methods": modes,
                "confidence": best_confidence,
                "evidence": evidence,
            }
        )
    return {
        "scope": VerdictScope.COMPLEX,
        "alternatives": alternatives,
        "ambiguity_reason": reason,
    }


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
        "scope": _scope_for_species(species),
        "alternatives": [],
        "ambiguity_reason": None,
    }


def _verdict_for(results: dict[str, Any]) -> dict[str, Any]:
    if len(results) > 1:
        return _arbitrate(results)
    first = next(iter(results.values()), {})
    if not first or "error" in first:
        return _verdict_payload("Unknown", "low", None, list(results), [], "0/1")
    species, conf, res = _species_of(first)
    return _verdict_payload(species, conf, res.get("strain"), list(results), [], "1/1")


def _identify_one(
    contigs: str | list[str],
    mode: str,
    db_dir: str | None,
    threads: int = 4,
    marker_db: str | None = None,
    db_selections: dict[str, list[Path]] | None = None,
) -> dict[str, Any]:
    files = list(contigs) if isinstance(contigs, (list, tuple)) else [contigs]
    primary = files[0]
    reads: str | list[str] = files if len(files) > 1 else primary
    results: dict[str, Any] = {}
    selections = db_selections or {}

    modes = _expand_modes(mode)
    unknown = [m for m in modes if m not in _ISOLATE_MODES]
    if unknown:
        raise ValueError(f"unsupported isolate mode(s): {', '.join(unknown)}")
    for m in modes:
        selected_paths = selections.get(m, [])
        try:
            if m == "marker":
                from gside.analysis.multigene_identifier import identify_multigene

                marker_paths: list[str] = [str(p) for p in selected_paths]
                if not marker_paths:
                    marker_path = marker_db
                    if marker_path is None and db_dir:
                        for candidate in (
                            Path(db_dir) / "D1_marker" / "default",
                            Path(db_dir) / "D1_marker",
                        ):
                            if candidate.is_dir():
                                marker_path = str(candidate)
                                break
                    if marker_path is not None:
                        marker_paths = [marker_path]
                if marker_paths:
                    results[m] = identify_multigene(
                        reads,
                        threads=threads,
                        marker_db=marker_paths[0] if len(marker_paths) == 1 else marker_paths,
                    ).to_dict()
                else:
                    results[m] = identify_multigene(reads, threads=threads).to_dict()
            elif m in ("ani", "mash"):
                from gside.analysis.ani_identifier import identify_by_ani

                results[m] = identify_by_ani(
                    reads,
                    mode=m,
                    db_dir=db_dir,
                    threads=threads,
                    dbs=selected_paths or None,
                ).to_dict()
            elif m == "sourmash":
                from gside.analysis.sourmash_identifier import (
                    identify_by_sourmash as identify_sourmash,
                )

                results[m] = identify_sourmash(
                    reads, db_dir=db_dir, dbs=selected_paths or None
                ).to_dict()
            elif m == "gtdbtk":
                from gside.analysis.gtdbtk_identifier import identify_by_gtdbtk

                results[m] = identify_by_gtdbtk(primary, threads=threads).to_dict()
        except Exception as e:  # noqa: BLE001 — CLI 边界：单法失败不终止其他方法
            entry = {"method": m, "error": str(e)[:300]}
            entry.update(_structured_missing_db_fields(entry["error"]))
            results[m] = entry

    verdict = _verdict_for(results)

    payload = {
        "analysis_type": "species_identification",
        "tool": "gside",
        "version": __version__,
        "contigs": primary,
        "methods": results,
        "verdict": verdict,
    }
    return payload


_LAYER = {"gtdbtk": 3, "ani": 2, "mash": 2, "sourmash": 2, "marker": 1}
_CONF_RANK = {"high": 2, "medium": 1, "low": 0}


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
    verdict = _verdict_payload(species, conf, res.get("strain"), [mode], supporting, agreement)
    verdict.update(_verdict_semantics(species, successful, results))
    return verdict


if __name__ == "__main__":
    main()
