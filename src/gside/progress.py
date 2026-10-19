"""Shared stderr-bound rich progress bar for multi-file runs (gmlst design).

Modeled on gmlst's ``make_progress``: spinner, description, bar, M-of-N,
percent, elapsed time — on a console bound to stderr so bar rendering can
never pollute stdout data streams (gside's stdout-purity law, same as
gapit). ``redirect_stdout=False`` / ``redirect_stderr=False`` are pinned
so rich's Live cannot reroute streamed report bytes.

The bar is ACTIVE only when ``--quiet`` is off, more than one input file
is given, and stderr is a TTY. When inactive, :func:`species_progress`
yields a no-op shim (task id ``None``) so call sites keep one code path.
"""

from __future__ import annotations

from collections.abc import Generator, Sequence
from contextlib import contextmanager
from typing import Final, Protocol

from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)

STATUS_CONSOLE: Final = Console(stderr=True)


class BatchProgress(Protocol):
    def describe(self, text: str) -> None: ...
    def advance(self) -> None: ...


class IdleBatch:
    def describe(self, text: str) -> None:
        return None

    def advance(self) -> None:
        return None


class _LiveBatch:
    def __init__(self, progress: Progress, task_id: TaskID) -> None:
        self._progress = progress
        self._task_id = task_id

    def describe(self, text: str) -> None:
        self._progress.update(self._task_id, description=text)

    def advance(self) -> None:
        self._progress.advance(self._task_id)


_IDLE: Final = IdleBatch()


def make_progress() -> Progress:
    """The gmlst column set on the stderr console; redirects pinned off."""
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=STATUS_CONSOLE,
        redirect_stdout=False,
        redirect_stderr=False,
    )


@contextmanager
def species_progress(
    inputs: Sequence[object], *, quiet: bool, description: str = "Identifying"
) -> Generator[tuple[BatchProgress, TaskID | None], None, None]:
    """Yield ``(bar, task_id)`` for a multi-file species run.

    ``task_id`` is None when the bar is INACTIVE — under ``--quiet``, when
    only one file is given, or when stderr is not a TTY — and then ``bar``
    is the no-op shim.
    """
    progress = None if quiet or len(inputs) <= 1 else make_progress()
    if progress is not None and progress.console.is_terminal:
        with progress:
            task_id = progress.add_task(description, total=len(inputs))
            yield _LiveBatch(progress, task_id), task_id
        return
    yield _IDLE, None
