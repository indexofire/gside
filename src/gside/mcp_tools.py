"""MCP tool implementations (protocol: gside.mcp; declarations: gside.mcp_schemas).

Every tool mirrors its CLI twin by calling the SAME shared callables —
never a reimplementation. Handlers validate their arguments and raise
ValueError with a clear message; the protocol layer renders failures as
isError envelopes. Parsed JSON-RPC argument containers stay dict[str, Any]
until the per-field isinstance checks below pin concrete types (the
sanctioned Any boundary, cli.py precedent).
"""

from __future__ import annotations

import contextlib
import json
import sys
from collections.abc import Callable
from typing import Any

_FORMATS = ("table", "tsv", "csv", "json")


def _contigs(arguments: dict[str, Any]) -> list[str]:
    value: Any = arguments.get("contigs", [])
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError("contigs must be an array of file path strings")
    if not value:
        raise ValueError("contigs must contain at least one file path")
    # value is proven list[str] by the checks above
    return list(value)


def _reads(arguments: dict[str, Any]) -> list[str]:
    value: Any = arguments.get("reads", [])
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError("reads must be an array of file path strings")
    if not value:
        raise ValueError("reads must contain at least one file path")
    if len(value) > 2:
        raise ValueError("reads accepts at most 2 files (1 = single-end/contigs, 2 = paired)")
    return list(value)


def _string(arguments: dict[str, Any], key: str, default: str) -> str:
    value = arguments.get(key, default)
    if not isinstance(value, str):
        raise ValueError(f"{key} must be a string")
    return value


def _optional_string(arguments: dict[str, Any], key: str) -> str | None:
    value = arguments.get(key)
    if value is not None and not isinstance(value, str):
        raise ValueError(f"{key} must be a string")
    return value


def _required_string(arguments: dict[str, Any], key: str) -> str:
    value = arguments.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{key} must be a non-empty string")
    return value


def _threads(arguments: dict[str, Any]) -> int:
    value = arguments.get("threads", 4)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("threads must be an integer >= 1")
    return value


def _format(arguments: dict[str, Any]) -> str:
    value = arguments.get("format", "json")
    if not isinstance(value, str) or value not in _FORMATS:
        raise ValueError(f"format must be one of {', '.join(_FORMATS)}: got {value}")
    return value


def _tool_species(arguments: dict[str, Any]) -> str:
    contigs = _contigs(arguments)
    mode = _string(arguments, "mode", "marker")
    from gside.cli import MODES

    if mode not in MODES:
        raise ValueError(f"mode must be one of {', '.join(MODES)}: got {mode}")
    db_dir = _optional_string(arguments, "db_dir")
    threads = _threads(arguments)
    from gside import cli

    samples = cli._group_read_samples(contigs)
    payloads = [
        cli._identify_one(sample[0] if len(sample) == 1 else sample, mode, db_dir, threads=threads)
        for sample in samples
    ]
    report: Any = payloads[0] if len(payloads) == 1 else payloads
    return json.dumps(report, ensure_ascii=False, indent=2)


def _tool_meta(arguments: dict[str, Any]) -> str:
    reads = _reads(arguments)
    mode = _string(arguments, "mode", "all")
    from gside.cli import META_MODES

    if mode not in META_MODES:
        raise ValueError(f"mode must be one of {', '.join(META_MODES)}: got {mode}")
    db_dir = _optional_string(arguments, "db_dir")
    threads = _threads(arguments)
    from gside import cli

    return json.dumps(cli._meta_payload(reads, mode, db_dir, threads), ensure_ascii=False, indent=2)


def _tool_db_status(arguments: dict[str, Any]) -> str:
    from gside.db import format_status

    return format_status(_format(arguments))


def _tool_db_list(arguments: dict[str, Any]) -> str:
    from gside.db import format_tiers

    return format_tiers(_format(arguments))


def _tool_db_setup(arguments: dict[str, Any]) -> str:
    tier = _required_string(arguments, "tier")
    from gside.cli import TIERS

    if tier not in TIERS:
        raise ValueError(f"tier must be one of {', '.join(TIERS)}: got {tier}")
    source = _optional_string(arguments, "source") or ""
    from gside.db import db_setup

    # db_setup prints progress; MCP stdout is protocol-only, so route it to stderr
    with contextlib.redirect_stdout(sys.stderr):
        results = db_setup(tier, source)
    return json.dumps(results, ensure_ascii=False, indent=2)


TOOL_HANDLERS: dict[str, Callable[[dict[str, Any]], str]] = {
    "species": _tool_species,
    "meta": _tool_meta,
    "db_status": _tool_db_status,
    "db_list": _tool_db_list,
    "db_setup": _tool_db_setup,
}
