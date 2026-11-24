"""Unit tests pinning explicit state machines across gside (UX/DX/AX).

Pure, hermetic, fast: no binaries, no network, no external databases.
"""

from __future__ import annotations

import hashlib
import io
import json
import threading
import time
from io import StringIO
from pathlib import Path
from typing import Any

import pytest

from gside import cli, db
from gside.analysis._thresholds import Confidence, ani_confidence
from gside.analysis.ani_identifier import _mash_confidence
from gside.analysis.multigene_identifier import (
    ControlSource,
    DepthGateMode,
    _apply_depth_gate,
)
from gside.analysis.taxonomic_validator import (
    QualityVerdict,
    ToolStatus,
    ValidationStatus,
    _aggregate_status,
    _quality_verdict,
)
from gside.cli import (
    InputKind,
    VerdictScope,
    _AutoSetupCoordinator,
    _input_kind,
    _json_stream,
    _OutputFileWriter,
    _run_species,
    _scope_for_species,
    _TierProvisionState,
    _WriterState,
)
from gside.db_registry import layout_valid
from gside.mcp import McpSession, McpState, _Frame, serve
from gside.mcp_schemas import TOOLS
from gside.progress import IdleBatch, ProgressState, _LiveBatch


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


class TestProgressStateMachine:
    def test_live_batch_nested_suspend_stops_and_starts_once(self) -> None:
        progress = _FakeProgress()
        bar = _LiveBatch(progress, 0)
        assert bar.state is ProgressState.RUNNING
        with bar.suspend():
            assert bar.state is ProgressState.SUSPENDED
            with bar.suspend():
                assert bar.state is ProgressState.SUSPENDED
            assert bar.state is ProgressState.SUSPENDED
            assert progress.events == ["stop"]
        assert bar.state is ProgressState.RUNNING
        assert progress.events == ["stop", "start"]

    def test_live_batch_closed_does_not_restart_progress(
        self, monkeypatch: Any, capsys: Any
    ) -> None:
        progress = _FakeProgress()
        bar = _LiveBatch(progress, 0)
        bar.close()
        assert bar.state is ProgressState.CLOSED

        buf = io.StringIO()
        monkeypatch.setattr("sys.stdout", buf)
        bar.write("after-close\n")
        bar.notice("notice-after-close")
        with bar.suspend():
            pass

        assert buf.getvalue() == "after-close\n"
        assert "notice-after-close" in capsys.readouterr().err
        assert progress.events == []

    def test_idle_batch_state_is_idle_and_suspend_is_noop(self) -> None:
        bar = IdleBatch()
        assert bar.state is ProgressState.IDLE
        with bar.suspend():
            assert bar.state is ProgressState.IDLE


class TestOutputFileWriterStateMachine:
    def test_writer_lifecycle_idle_open_closed(self, tmp_path: Path) -> None:
        target = tmp_path / "out.tsv"
        writer = _OutputFileWriter(target)
        assert writer.state is _WriterState.IDLE
        assert not target.exists()

        writer("a\n")
        assert writer.state is _WriterState.OPEN
        assert target.read_text() == "a\n"

        writer.close()
        assert writer.state is _WriterState.CLOSED
        # Closing twice is idempotent
        writer.close()
        assert writer.state is _WriterState.CLOSED

    def test_writer_rejects_write_after_close(self, tmp_path: Path) -> None:
        writer = _OutputFileWriter(tmp_path / "out.tsv")
        writer.close()
        with pytest.raises(RuntimeError, match="closed"):
            writer("late\n")

    def test_writer_context_manager_closes_on_exit(self, tmp_path: Path) -> None:
        target = tmp_path / "ctx.tsv"
        with _OutputFileWriter(target) as writer:
            writer("row\n")
            assert writer.state is _WriterState.OPEN
        assert writer.state is _WriterState.CLOSED
        assert target.read_text() == "row\n"

    def test_json_stream_empty_single_yields_nothing(self) -> None:
        assert list(_json_stream(iter([]), single=True)) == []


def _missing_db_payload(
    contigs: str, method: str = "sourmash", tier: str = "sourmash"
) -> dict[str, Any]:
    return {
        "analysis_type": "species_identification",
        "tool": "gside",
        "contigs": contigs,
        "methods": {
            method: {
                "method": method,
                "error": f"{method} DB not found (run: gside db setup --tier {tier})",
                "suggested_command": f"gside db setup --tier {tier}",
                "missing_tier": tier,
            }
        },
        "verdict": {"species": "Unknown", "confidence": "low", "basis": [method]},
    }


def _ok_payload(
    contigs: str, method: str = "sourmash", species: str = "Escherichia_coli"
) -> dict[str, Any]:
    return {
        "analysis_type": "species_identification",
        "tool": "gside",
        "contigs": contigs,
        "methods": {
            method: {
                "method": method,
                "result": {"species": species, "confidence": "high"},
            }
        },
        "verdict": {"species": species, "confidence": "high", "basis": [method]},
    }


class TestAutoSetupCoordinatorStateMachine:
    def test_batch_deduplicates_tier_setup_across_samples(
        self, monkeypatch: Any, capsys: Any
    ) -> None:
        setups: list[str] = []
        installed = False

        def fake_identify(
            contigs: Any,
            mode: str,
            db_dir: Any,
            threads: int = 4,
            marker_db: Any = None,
            db_selections: Any = None,
        ) -> dict[str, Any]:
            if not installed:
                return _missing_db_payload(str(contigs))
            return _ok_payload(str(contigs))

        def fake_setup(tier: str) -> None:
            nonlocal installed
            setups.append(tier)
            installed = True

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", fake_setup)

        assert _run_species(["a.fna", "b.fna", "c.fna"], "sourmash", None, auto_setup=True) == 0
        payloads = json.loads(capsys.readouterr().out)
        assert setups == ["sourmash"]
        assert [p["verdict"]["species"] for p in payloads] == [
            "Escherichia_coli",
            "Escherichia_coli",
            "Escherichia_coli",
        ]

    def test_failed_tier_does_not_retry_for_subsequent_samples(
        self, monkeypatch: Any, capsys: Any
    ) -> None:
        setups: list[str] = []

        def fake_identify(
            contigs: Any,
            mode: str,
            db_dir: Any,
            threads: int = 4,
            marker_db: Any = None,
            db_selections: Any = None,
        ) -> dict[str, Any]:
            return _missing_db_payload(str(contigs))

        def boom(tier: str) -> None:
            setups.append(tier)
            raise RuntimeError("offline")

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", boom)

        assert _run_species(["a.fna", "b.fna", "c.fna"], "sourmash", None, auto_setup=True) == 0
        payloads = json.loads(capsys.readouterr().out)
        assert setups == ["sourmash"]
        assert all(p["verdict"]["species"] == "Unknown" for p in payloads)

    def test_setup_error_dict_marks_tier_failed_without_retry(
        self, monkeypatch: Any, capsys: Any
    ) -> None:
        import gside.db as dbmod

        calls: list[str] = []

        def fake_identify(
            contigs: Any,
            mode: str,
            db_dir: Any,
            threads: int = 4,
            marker_db: Any = None,
            db_selections: Any = None,
        ) -> dict[str, Any]:
            calls.append(f"{contigs}:{mode}")
            return _missing_db_payload(str(contigs))

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(
            dbmod, "db_setup", lambda tier, source="": {"D4_sourmash": "ERROR: down"}
        )

        assert _run_species(["a.fna", "b.fna"], "sourmash", None, auto_setup=True) == 0
        json.loads(capsys.readouterr().out)
        # No futile retry for a.fna or b.fna after db_setup returned ERROR
        assert calls == ["a.fna:sourmash", "b.fna:sourmash"]

    def test_parallel_jobs_serialize_single_tier_install(
        self, monkeypatch: Any, capsys: Any
    ) -> None:
        setup_count = 0
        setup_lock = threading.Lock()
        installed = threading.Event()
        started = threading.Barrier(3)

        def fake_identify(
            contigs: Any,
            mode: str,
            db_dir: Any,
            threads: int = 4,
            marker_db: Any = None,
            db_selections: Any = None,
        ) -> dict[str, Any]:
            if not installed.is_set():
                started.wait(timeout=5)
                return _missing_db_payload(str(contigs))
            return _ok_payload(str(contigs))

        def slow_setup(tier: str) -> None:
            nonlocal setup_count
            time.sleep(0.05)
            with setup_lock:
                setup_count += 1
            installed.set()

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", slow_setup)

        assert (
            _run_species(["a.fna", "b.fna", "c.fna"], "sourmash", None, jobs=3, auto_setup=True)
            == 0
        )
        payloads = json.loads(capsys.readouterr().out)
        assert setup_count == 1
        assert [p["verdict"]["species"] for p in payloads] == [
            "Escherichia_coli",
            "Escherichia_coli",
            "Escherichia_coli",
        ]

    def test_coordinator_suspends_live_progress_bar(self, monkeypatch: Any) -> None:
        progress = _FakeProgress()
        bar = _LiveBatch(progress, 0)
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: None)

        coordinator = _AutoSetupCoordinator()
        assert coordinator.state_of("sourmash") is _TierProvisionState.PENDING
        assert coordinator.provision("sourmash", bar=bar) is True
        assert coordinator.state_of("sourmash") is _TierProvisionState.INSTALLED
        assert progress.events == ["stop", "start"]

    def test_auto_setup_retry_preserves_db_selections(self, monkeypatch: Any) -> None:
        seen_selections: list[Any] = []

        def fake_identify(
            contigs: Any,
            mode: str,
            db_dir: Any,
            threads: int = 4,
            marker_db: Any = None,
            db_selections: Any = None,
        ) -> dict[str, Any]:
            seen_selections.append(db_selections)
            return _ok_payload(str(contigs))

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: None)

        selections = {"sourmash": [Path("/tmp/custom")]}
        payload = _missing_db_payload("a.fna")
        cli._auto_setup_retry("a.fna", payload, "sourmash", None, 4, db_selections=selections)
        assert seen_selections == [selections]

    def test_auto_setup_retry_preserves_marker_db(self, monkeypatch: Any) -> None:
        seen_marker_dbs: list[Any] = []

        def fake_identify(
            contigs: Any,
            mode: str,
            db_dir: Any,
            threads: int = 4,
            marker_db: Any = None,
        ) -> dict[str, Any]:
            seen_marker_dbs.append(marker_db)
            return _ok_payload(str(contigs))

        monkeypatch.setattr(cli, "_identify_one", fake_identify)
        monkeypatch.setattr(cli, "_setup_missing_tier", lambda tier: None)

        payload = _missing_db_payload("a.fna")
        cli._auto_setup_retry("a.fna", payload, "sourmash", None, 4, marker_db="/tmp/markers")
        assert seen_marker_dbs == ["/tmp/markers"]


class TestMcpSessionStateMachine:
    def test_lifecycle_uninitialized_initialized_ready_closed(self) -> None:
        session = McpSession()
        assert session.state is McpState.UNINITIALIZED

        resp = session.handle(
            _Frame.model_validate({"id": 1, "method": "initialize", "params": {}})
        )
        assert resp is not None and "result" in resp
        assert session.state is McpState.INITIALIZED

        notif = session.handle(_Frame.model_validate({"method": "notifications/initialized"}))
        assert notif is None
        assert session.state is McpState.READY

        session.close()
        assert session.state is McpState.CLOSED

    def test_ping_returns_empty_result_over_the_wire(self) -> None:
        raw = (
            json.dumps({"jsonrpc": "2.0", "id": 1, "method": "ping"})
            + "\n"
            + json.dumps({"jsonrpc": "2.0", "id": 2, "method": "initialize", "params": {}})
            + "\n"
            + json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"})
            + "\n"
            + json.dumps({"jsonrpc": "2.0", "id": 3, "method": "ping"})
            + "\n"
        )
        out = StringIO()
        serve(StringIO(raw), out)
        responses = [json.loads(line) for line in out.getvalue().splitlines()]
        assert responses[0] == {"jsonrpc": "2.0", "id": 1, "result": {}}
        assert responses[1]["id"] == 2
        assert responses[2] == {"jsonrpc": "2.0", "id": 3, "result": {}}

    def test_species_tool_description_is_well_formed(self) -> None:
        (species,) = [tool for tool in TOOLS if tool["name"] == "species"]
        desc = str(species["description"])
        assert "one methods" not in desc
        assert "FASTA" in desc or "assemblies" in desc


class TestDomainStateEnums:
    def test_taxonomic_validator_enums_and_aggregate_status(self) -> None:
        assert _aggregate_status([ToolStatus.OK, ToolStatus.OK]) is ValidationStatus.OK
        assert _aggregate_status([ToolStatus.OK, ToolStatus.MISSING]) is ValidationStatus.PARTIAL
        assert _aggregate_status([ToolStatus.MISSING, ToolStatus.FAILED]) is ValidationStatus.ERROR
        assert _quality_verdict(99.0, 1.0) is QualityVerdict.PASS
        assert _quality_verdict(80.0, 1.0) is QualityVerdict.CHECK
        assert _quality_verdict(40.0, 1.0) is QualityVerdict.WARNING
        assert _quality_verdict(None, 1.0) is QualityVerdict.UNKNOWN

    def test_confidence_enum_and_deduplicated_mash_confidence(self) -> None:
        assert ani_confidence(96.0) is Confidence.HIGH
        assert ani_confidence(94.0) is Confidence.MEDIUM
        assert ani_confidence(90.0) is Confidence.LOW
        assert _mash_confidence(0.98) is Confidence.HIGH
        assert _mash_confidence(0.92) is Confidence.MEDIUM
        assert _mash_confidence(0.85) is Confidence.LOW

    def test_multigene_depth_gate_mode_transitions(self) -> None:
        # Assembly hits (no mean_depth) -> DISABLED
        asm_gate = _apply_depth_gate(
            [{"gene": "inva", "identity": 99.0, "coverage": 100.0}], [], 5.0, 0.10
        )
        assert asm_gate.mode is DepthGateMode.DISABLED
        assert asm_gate.control_source is ControlSource.NONE

        # Reads hit with no passing control -> ABSOLUTE_ONLY
        abs_gate = _apply_depth_gate(
            [{"gene": "inva", "identity": 80.0, "coverage": 50.0, "mean_depth": 10.0}],
            [],
            5.0,
            0.10,
        )
        assert abs_gate.mode is DepthGateMode.ABSOLUTE_ONLY
        assert abs_gate.control_source is ControlSource.NONE

        # Reads hit with dedicated control -> RELATIVE + DEDICATED
        rel_gate = _apply_depth_gate(
            [
                {"gene": "rpob", "identity": 99.0, "coverage": 100.0, "mean_depth": 50.0},
                {"gene": "inva", "identity": 99.0, "coverage": 100.0, "mean_depth": 40.0},
            ],
            ["rpob"],
            5.0,
            0.10,
        )
        assert rel_gate.mode is DepthGateMode.RELATIVE
        assert rel_gate.control_source is ControlSource.DEDICATED

    def test_cli_input_kind_and_verdict_scope_enums(self) -> None:
        assert _input_kind("a.fna.gz") is InputKind.FASTA
        assert _input_kind("r1.fastq.gz") is InputKind.FASTQ
        assert _input_kind("readme.txt") is InputKind.UNKNOWN
        assert _scope_for_species("Escherichia_coli") is VerdictScope.SPECIES
        assert _scope_for_species("Salmonella") is VerdictScope.GENUS
        assert _scope_for_species("Salmonella_enterica") is VerdictScope.SPECIES
        assert _scope_for_species("Bacillus_cereus_group") is VerdictScope.COMPLEX
        assert _scope_for_species("Unknown") is VerdictScope.UNKNOWN


class TestDbDownloadAndLayoutStateMachine:
    def test_download_file_atomic_part_cleaned_on_failure(
        self, tmp_path: Path, monkeypatch: Any
    ) -> None:
        class _BrokenResponse:
            headers = {"Content-Length": "10"}

            def __init__(self) -> None:
                self._calls = 0

            def read(self, _size: int = -1) -> bytes:
                self._calls += 1
                if self._calls == 1:
                    return b"part"
                raise OSError("connection reset")

            def __enter__(self) -> _BrokenResponse:
                return self

            def __exit__(self, *args: Any) -> bool:
                return False

        monkeypatch.setattr(db.urllib.request, "urlopen", lambda *a, **k: _BrokenResponse())
        target = tmp_path / "mash.msh.gz"
        with pytest.raises(OSError, match="connection reset"):
            db._download_file("https://example.invalid/mash.msh.gz", target)

        assert not target.exists()
        assert not target.with_name(f"{target.name}.part").exists()

    def test_sourmash_lineages_mismatch_cleans_both_artifacts(
        self, tmp_path: Path, monkeypatch: Any
    ) -> None:
        root = tmp_path / "db"
        monkeypatch.setattr(db, "SPECIES_DB_DIR", root)
        monkeypatch.setattr(
            db, "_download_file", lambda url, dst: Path(dst).write_bytes(b"payload")
        )
        monkeypatch.setattr(db, "SOURMASH_SIG_SHA256", hashlib.sha256(b"payload").hexdigest())
        monkeypatch.setattr(db, "SOURMASH_LINEAGES_SHA256", "0" * 64)

        dst = root / "D4_sourmash"
        dst.mkdir(parents=True)
        msg = db._download_sourmash_farm(dst)
        assert msg.startswith("ERROR: SHA256 mismatch")
        # Both staging files are cleaned up on verification failure
        assert not (dst / "default" / "gtdb-rs226-reps.k31.sig.zip").exists()
        assert not (dst / "default" / "gtdb-rs226-reps.lineages.csv").exists()

    def test_db_registry_accepts_mash_payload_bin(self, tmp_path: Path) -> None:
        (tmp_path / "payload.bin").write_bytes(b"msh")
        assert layout_valid("mash", tmp_path) is True
