"""Single home for binary discovery + checked subprocess execution (engine backends)."""

from __future__ import annotations

import subprocess
from collections.abc import Callable

from ..config import pixi_path, which

__all__ = ["pixi_path", "require_bin", "run_checked", "which"]


def require_bin(
    tool: str,
    *,
    hint: str = "",
    resolver: Callable[[str], str | None] | None = None,
) -> str:
    """Resolve ``tool`` to a path or raise RuntimeError with ``hint`` as the message.

    Backends pass their module-level ``which`` as ``resolver`` so per-backend
    monkeypatching of that name keeps working.
    """
    found = (resolver or which)(tool)
    if not found:
        raise RuntimeError(hint or f"{tool} not found in PATH")
    return found


def run_checked(
    cmd: list[str],
    *,
    timeout: int,
    name: str,
    exit_in_msg: bool = True,
    truncate: int | None = 500,
) -> subprocess.CompletedProcess[str]:
    """Run ``cmd`` captured; on nonzero exit raise RuntimeError as
    ``"{name}{' (exit N)' if exit_in_msg}: {stderr.strip()[:truncate]}"``.
    """
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        detail = proc.stderr.strip()
        if truncate is not None:
            detail = detail[:truncate]
        exit_part = f" (exit {proc.returncode})" if exit_in_msg else ""
        raise RuntimeError(f"{name}{exit_part}: {detail}")
    return proc
