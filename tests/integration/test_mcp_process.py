"""Subprocess regression: `python -m gside.mcp` must boot standalone.

The unit suite imports gside.cli before gside.mcp, which masks import-order
cycles; this test enters through the module itself, exactly like an MCP
client spawning `gside-mcp` / `python -m gside.mcp`.
"""

from __future__ import annotations

import json
import subprocess
import sys


def test_module_main_answers_initialize() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "gside.mcp"],
        input='{"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {}}\n',
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    )
    (response,) = [json.loads(line) for line in proc.stdout.splitlines()]
    assert response["id"] == 1
    assert response["result"]["serverInfo"]["name"] == "gside"
    assert response["result"]["protocolVersion"] == "2025-06-18"
