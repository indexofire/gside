"""Progress-bar/output coordination: rows must not be erased by the Live bar."""

from __future__ import annotations

import io
from typing import Any

import pytest

from gside import cli
from gside.progress import IdleBatch, _LiveBatch, make_progress


class _FakeProgress:
    def __init__(self) -> None:
        self.events: list[str] = []

    def stop(self) -> None:
        self.events.append("stop")

    def start(self) -> None:
        self.events.append("start")

    def update(self, task_id: Any, **kwargs: Any) -> None:
        self.events.append(f"update:{task_id}:{kwargs.get('description', '')}")

    def advance(self, task_id: Any) -> None:
        self.events.append(f"advance:{task_id}")


class _RecordingBar:
    def __init__(self) -> None:
        self.writes: list[str] = []
        self.notices: list[str] = []

    def describe(self, text: str) -> None:
        pass

    def advance(self) -> None:
        pass

    def write(self, text: str) -> None:
        self.writes.append(text)

    def notice(self, text: str) -> None:
        self.notices.append(text)


class TestLiveBatch:
    def test_describe_and_advance_forward_to_progress(self):
        progress = _FakeProgress()
        bar = _LiveBatch(progress, 7)
        bar.describe("Identifying (2 jobs)")
        bar.advance()
        assert progress.events == ["update:7:Identifying (2 jobs)", "advance:7"]

    def test_write_suspends_and_restarts_the_bar(self, monkeypatch):
        progress = _FakeProgress()
        buffer = io.StringIO()
        monkeypatch.setattr("sys.stdout", buffer)
        _LiveBatch(progress, 0).write("row\n")
        assert buffer.getvalue() == "row\n"
        assert progress.events == ["stop", "start"]

    def test_notice_suspends_and_restarts_the_bar(self, monkeypatch, capsys):
        progress = _FakeProgress()
        _LiveBatch(progress, 0).notice("Results written to out.tsv")
        assert progress.events == ["stop", "start"]
        assert "Results written" in capsys.readouterr().err

    def test_start_is_restored_when_the_write_fails(self, monkeypatch):
        progress = _FakeProgress()

        class Broken:
            def write(self, text: str) -> int:
                raise OSError("pipe closed")

            def flush(self) -> None:
                pass

        monkeypatch.setattr("sys.stdout", Broken())
        with pytest.raises(OSError):
            _LiveBatch(progress, 0).write("row\n")
        assert progress.events == ["stop", "start"]


class TestIdleBatch:
    def test_write_goes_to_stdout(self, capsys):
        IdleBatch().write("row\n")
        assert capsys.readouterr().out == "row\n"

    def test_notice_goes_to_stderr(self, capsys):
        IdleBatch().notice("notice")
        assert capsys.readouterr().err == "notice\n"


class TestEmissionUsesBar:
    def test_streamed_chunks_go_through_bar(self):
        bar = _RecordingBar()
        assert cli._emit_chunks(iter(["a\n", "b\n"]), None, bar) == 0
        assert bar.writes == ["a\n", "b\n"]

    def test_json_document_goes_through_bar(self):
        bar = _RecordingBar()
        assert cli._emit_species(iter([{"x": 1}]), 1, "json", None, bar=bar) == 0
        assert bar.writes and '"x": 1' in bar.writes[0]

    def test_file_output_notice_goes_through_bar(self, tmp_path):
        bar = _RecordingBar()
        target = tmp_path / "out.tsv"
        assert cli._emit_chunks(iter(["a\n"]), str(target), bar) == 0
        assert target.read_text() == "a\n"
        assert bar.writes == []
        assert bar.notices and "Results written" in bar.notices[0]


def test_progress_is_transient_so_stop_clears_the_bar():
    progress = make_progress()
    assert progress.live.transient is True
    progress.stop()
