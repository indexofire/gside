r"""Minimal MCP stdio server: gside as a tool provider for agents.

Hand-rolled on purpose (AGENTS.md keeps the dependency list deliberately
short): no `mcp` SDK — just the essential MCP stdio behavior, JSON-RPC 2.0,
one message per line on stdin and one response line on stdout. Limitations:
single messages only (no batch arrays); non-JSON lines are ignored silently
(robustness over -32700); JSON lines that fail frame validation get a
-32600/-32602 error when they carry a request id, so no client ever hangs.
This module is the PROTOCOL only — frames, dispatch, and the serve loop;
the tool implementations live in :mod:`gside.mcp_tools` and the tools/list
declarations (inputSchemas) in :mod:`gside.mcp_schemas` (species, db_status,
db_list, db_setup). Tool failures return isError=true with the exception
rendered as text; stdout is protocol-only.
"""

# No `from __future__ import annotations` here on purpose: pydantic resolves
# BaseModel field annotations through the defining module's namespace, and
# this module also runs as __main__ (`python -m gside.mcp`), where postponed
# string annotations cannot be resolved.
#
# TOOLS/TOOL_HANDLERS are imported lazily (at dispatch time), not here:
# gside.cli imports this module when registering `gside mcp`, and a
# module-level import would close the cycle cli -> mcp -> mcp_schemas -> cli
# and break `python -m gside.mcp` as the first entry point.
import json
import sys
from collections.abc import Iterable
from typing import Any, TextIO, TypeGuard

import click
from pydantic import BaseModel, ConfigDict, ValidationError

from gside import __version__

# Parsed JSON-RPC values are the one sanctioned Any boundary (cli.py
# precedent): frame/argument containers stay dict[str, Any] until the
# per-field checks in mcp_tools pin concrete types.
JsonRpcId = str | int | float | None

PROTOCOL_VERSION = "2025-06-18"


class _Frame(BaseModel, frozen=True):
    """One stdin frame (request or notification); ``id`` presence in
    model_fields_set distinguishes requests, per JSON-RPC 2.0."""

    model_config = ConfigDict(extra="ignore")

    method: str | None = None
    params: dict[str, Any] | None = None
    id: JsonRpcId = None


class _ToolCall(BaseModel, frozen=True):
    """tools/call params: the tool ``name`` and its ``arguments`` object."""

    model_config = ConfigDict(extra="ignore")

    name: str = ""
    arguments: dict[str, Any] | None = None


def _result(id_value: JsonRpcId, result: dict[str, object]) -> dict[str, object]:
    return {"jsonrpc": "2.0", "id": id_value, "result": result}


def _error(id_value: JsonRpcId, code: int, message: str) -> dict[str, object]:
    return {"jsonrpc": "2.0", "id": id_value, "error": {"code": code, "message": message}}


def _text(text: str, *, is_error: bool) -> dict[str, object]:
    return {"content": [{"type": "text", "text": text}], "isError": is_error}


def _initialize(params: dict[str, Any] | None) -> dict[str, object]:
    version = (params or {}).get("protocolVersion", "")
    requested = version if isinstance(version, str) else ""
    return {
        "protocolVersion": requested or PROTOCOL_VERSION,
        "capabilities": {"tools": {}},
        "serverInfo": {"name": "gside", "version": __version__},
    }


def _tools_call(id_value: JsonRpcId, params: dict[str, Any] | None) -> dict[str, object]:
    from gside.mcp_tools import TOOL_HANDLERS

    try:
        call = _ToolCall.model_validate(params or {})
    except ValidationError:
        return _error(id_value, -32602, "invalid tools/call params")
    tool = TOOL_HANDLERS.get(call.name)
    if tool is None:
        return _error(id_value, -32602, f"unknown tool: {call.name or '(none)'}")
    try:
        text = tool(call.arguments or {})
    except Exception as exc:  # cli.py boundary semantics: envelope, not a crash
        return _result(id_value, _text(f"{type(exc).__name__}: {exc}", is_error=True))
    return _result(id_value, _text(text, is_error=False))


def _is_json_object(value: Any) -> TypeGuard[dict[str, Any]]:
    """Narrow a parsed JSON value to the sanctioned dict[str, Any] boundary
    (plain isinstance would surface dict[Unknown, Unknown])."""
    return isinstance(value, dict)


def _reject(value: Any) -> dict[str, object] | None:
    """A line that parsed as JSON but failed _Frame validation: answer the
    request id with -32600 (missing/non-string method) or -32602 (params not
    an object); id-less messages stay ignored like notifications."""
    if not _is_json_object(value):
        return None
    id_value = value.get("id")
    if id_value is None:
        return None
    if not isinstance(value.get("method"), str):
        return _error(id_value, -32600, "invalid request: method must be a string")
    return _error(id_value, -32602, "invalid request: params must be an object")


def _handle(frame: _Frame) -> dict[str, object] | None:
    """One parsed frame -> one response object; None = no response."""
    if "id" not in frame.model_fields_set:
        return None  # notifications never get responses
    if frame.method is None:
        return _error(frame.id, -32600, "invalid request: method is required")
    id_value = frame.id
    if frame.method == "initialize":
        return _result(id_value, _initialize(frame.params))
    if frame.method == "tools/list":
        from gside.mcp_schemas import TOOLS

        return _result(id_value, {"tools": list(TOOLS)})
    if frame.method == "tools/call":
        return _tools_call(id_value, frame.params)
    return _error(id_value, -32601, f"method not found: {frame.method}")


def serve(stdin: Iterable[str], stdout: TextIO) -> None:
    """Serve newline-delimited JSON-RPC 2.0 until EOF; non-JSON lines and
    id-less invalid messages are ignored, requests that fail frame
    validation get a -32600/-32602 error instead of silence."""
    for line in stdin:
        try:
            value: object = json.loads(line)
        except json.JSONDecodeError:
            continue
        try:
            frame = _Frame.model_validate(value)
        except ValidationError:
            response = _reject(value)
        else:
            response = _handle(frame)
        if response is not None:
            stdout.write(json.dumps(response, separators=(",", ":")) + "\n")
            stdout.flush()


def register_mcp_command(app: click.Group) -> None:
    """Attach the `mcp` subcommand (cli.py stays import + registration only)."""
    from gside.cli import HELP_SETTINGS  # lazy: gside.cli imports this module first

    app.command("mcp", context_settings=HELP_SETTINGS)(main)


def main() -> None:
    """MCP stdio server: gside as a tool provider for agents."""
    serve(sys.stdin, sys.stdout)


if __name__ == "__main__":
    main()
