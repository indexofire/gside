"""Pure unit tests for the MCP stdio server (protocol, schemas, handlers).

No binaries, no databases, no network: species calls stub cli._identify_one,
db calls monkeypatch gside.db callables. Protocol tests feed newline-delimited
JSON-RPC lines to serve() and parse the response lines (gapit test recipe).
"""

from __future__ import annotations

import json
import runpy
from io import StringIO
from typing import Any

import pytest
from click.testing import CliRunner

from gside import __version__, cli
from gside.cli import MODES, TIERS, main
from gside.mcp import serve
from gside.mcp_schemas import TOOLS
from gside.mcp_tools import TOOL_HANDLERS


def exchange(*lines_or_messages: object) -> list[dict[str, Any]]:
    """Feed the loop one line per item (str = raw line, dict = JSON-RPC
    message); return the parsed response objects."""
    raw = "".join(
        (item if isinstance(item, str) else json.dumps(item)) + "\n" for item in lines_or_messages
    )
    out = StringIO()
    serve(StringIO(raw), out)
    return [json.loads(line) for line in out.getvalue().splitlines()]


def _payload(
    contigs: str = "a.fna", species: str = "Salmonella", conf: str = "high"
) -> dict[str, Any]:
    return {
        "analysis_type": "species_identification",
        "tool": "gside",
        "contigs": contigs,
        "methods": {"marker": {"species": species, "confidence": conf}},
        "verdict": {"species": species, "confidence": conf, "basis": ["marker"]},
    }


def _call(name: str, arguments: dict[str, Any]) -> dict[str, Any]:
    return {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "tools/call",
        "params": {"name": name, "arguments": arguments},
    }


_STATUS: dict[str, dict[str, Any]] = {
    "D1_marker": {"ready": True, "path": "/repo/data/db/D1_marker", "tier": "mini"},
    "D2_ani": {"ready": False, "path": "/repo/data/db/D2_ani", "tier": "panel"},
    "D3_mash": {"ready": False, "path": "/repo/data/db/D3_mash", "tier": "mash"},
    "D4_sourmash": {"ready": False, "path": "/repo/data/db/D4_sourmash", "tier": "sourmash"},
    "D5_kraken2": {"ready": False, "path": "/repo/data/db/D5_kraken2", "tier": "meta"},
    "D6_sylph": {"ready": False, "path": "/repo/data/db/D6_sylph", "tier": "meta_sylph"},
}


class TestServeProtocol:
    def test_initialize_advertises_gside_server_info(self):
        (response,) = exchange({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {}})
        result = response["result"]
        assert result["serverInfo"] == {"name": "gside", "version": __version__}
        assert result["capabilities"] == {"tools": {}}
        assert result["protocolVersion"] == "2025-06-18"

    def test_initialize_echoes_requested_protocol_version(self):
        (response,) = exchange(
            {
                "jsonrpc": "2.0",
                "id": 1,
                "method": "initialize",
                "params": {"protocolVersion": "2025-06-18"},
            }
        )
        assert response["result"]["protocolVersion"] == "2025-06-18"

    def test_initialize_non_string_version_falls_back(self):
        (response,) = exchange(
            {"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": 42}}
        )
        assert response["result"]["protocolVersion"] == "2025-06-18"

    def test_tools_list_advertises_every_tool(self):
        (response,) = exchange({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
        names = [tool["name"] for tool in response["result"]["tools"]]
        assert names == ["species", "meta", "db_status", "db_list", "db_setup"]

    def test_unknown_tool_gets_32602(self):
        (response,) = exchange(_call("nope", {}))
        assert response["id"] == 1
        assert response["error"]["code"] == -32602
        assert "unknown tool: nope" in response["error"]["message"]

    def test_missing_tool_name_gets_32602(self):
        (response,) = exchange({"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": {}})
        assert response["error"]["code"] == -32602
        assert "(none)" in response["error"]["message"]

    def test_non_string_tool_name_gets_32602(self):
        (response,) = exchange(
            {"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": {"name": 42}}
        )
        assert response["error"]["code"] == -32602
        assert response["error"]["message"] == "invalid tools/call params"

    def test_array_params_gets_32602_not_silence(self):
        (response,) = exchange(
            {"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": ["species"]}
        )
        assert response["id"] == 3
        assert response["error"]["code"] == -32602

    def test_unknown_method_gets_32601(self):
        (response,) = exchange({"jsonrpc": "2.0", "id": 5, "method": "resources/list"})
        assert response["error"]["code"] == -32601
        assert "resources/list" in response["error"]["message"]

    def test_request_without_method_gets_32600(self):
        (response,) = exchange({"jsonrpc": "2.0", "id": 1})
        assert response["id"] == 1
        assert response["error"]["code"] == -32600

    def test_non_string_method_gets_32600(self):
        (response,) = exchange({"jsonrpc": "2.0", "id": 4, "method": 42})
        assert response["id"] == 4
        assert response["error"]["code"] == -32600

    def test_non_json_line_is_ignored(self):
        responses = exchange("not json {", {"jsonrpc": "2.0", "id": 7, "method": "tools/list"})
        assert [r["id"] for r in responses] == [7]

    def test_bare_json_scalar_is_ignored(self):
        assert exchange("42") == []

    def test_notification_without_id_gets_no_response(self):
        assert exchange({"jsonrpc": "2.0", "method": "tools/list"}) == []

    def test_invalid_idless_messages_stay_ignored(self):
        responses = exchange(
            {"method": 42},
            {"jsonrpc": "2.0", "method": "tools/call", "params": ["species"]},
            {"jsonrpc": "2.0", "id": 9, "method": "tools/list"},
        )
        assert [r["id"] for r in responses] == [9]


class TestToolHandlers:
    def test_species_single_contig_returns_payload_dict(self, monkeypatch):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        payload = json.loads(TOOL_HANDLERS["species"]({"contigs": ["a.fna"]}))
        assert isinstance(payload, dict)
        assert payload["contigs"] == "a.fna"
        assert payload["verdict"]["species"] == "Salmonella"

    def test_species_multiple_contigs_returns_payload_list(self, monkeypatch):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        payloads = json.loads(TOOL_HANDLERS["species"]({"contigs": ["a.fna", "b.fna"]}))
        assert [p["contigs"] for p in payloads] == ["a.fna", "b.fna"]

    def test_species_forwards_mode_db_dir_threads(self, monkeypatch):
        seen: list[tuple[str, str, str | None, int]] = []

        def fake(contigs: str, mode: str, db_dir: str | None, threads: int = 4) -> dict[str, Any]:
            seen.append((contigs, mode, db_dir, threads))
            return _payload(contigs)

        monkeypatch.setattr(cli, "_identify_one", fake)
        TOOL_HANDLERS["species"](
            {"contigs": ["x.fna"], "mode": "panel", "db_dir": "custom-db", "threads": 8}
        )
        assert seen == [("x.fna", "panel", "custom-db", 8)]

    def test_db_status_returns_valid_json(self, monkeypatch):
        from gside import db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: _STATUS)
        payload = json.loads(TOOL_HANDLERS["db_status"]({}))
        assert payload == _STATUS

    def test_db_list_returns_valid_json(self, monkeypatch):
        from gside import db as dbmod

        monkeypatch.setattr(dbmod, "db_status", lambda: _STATUS)
        payload = json.loads(TOOL_HANDLERS["db_list"]({}))
        assert set(payload) == set(TIERS)
        assert payload["mini"]["ready"] is True
        assert payload["panel"]["ready"] is False

    def test_db_setup_returns_results_json(self, monkeypatch):
        from gside import db as dbmod

        monkeypatch.setattr(dbmod, "db_setup", lambda tier, source="": {"D2_ani": "already ready"})
        assert json.loads(TOOL_HANDLERS["db_setup"]({"tier": "panel"})) == {
            "D2_ani": "already ready"
        }

    def test_db_setup_error_result_still_json(self, monkeypatch):
        from gside import db as dbmod

        monkeypatch.setattr(dbmod, "db_setup", lambda tier, source="": {"D2_ani": "ERROR: boom"})
        assert json.loads(TOOL_HANDLERS["db_setup"]({"tier": "panel"})) == {"D2_ani": "ERROR: boom"}

    def test_db_setup_routes_progress_to_stderr(self, monkeypatch, capsys):
        from gside import db as dbmod

        def loud_setup(tier: str, source: str = "") -> dict[str, str]:
            print("progress line")
            return {"D2_ani": "ok"}

        monkeypatch.setattr(dbmod, "db_setup", loud_setup)
        TOOL_HANDLERS["db_setup"]({"tier": "panel"})
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "progress line" in captured.err

    def test_species_rejects_missing_contigs(self):
        with pytest.raises(ValueError, match="contigs"):
            TOOL_HANDLERS["species"]({})

    def test_species_rejects_non_list_contigs(self):
        with pytest.raises(ValueError, match="contigs must be an array"):
            TOOL_HANDLERS["species"]({"contigs": "a.fna"})

    def test_species_rejects_non_string_contig_items(self):
        with pytest.raises(ValueError, match="contigs must be an array"):
            TOOL_HANDLERS["species"]({"contigs": ["a.fna", 5]})

    def test_species_rejects_empty_contigs(self):
        with pytest.raises(ValueError, match="at least one"):
            TOOL_HANDLERS["species"]({"contigs": []})

    def test_species_rejects_unknown_mode(self):
        with pytest.raises(ValueError, match="mode must be one of"):
            TOOL_HANDLERS["species"]({"contigs": ["a.fna"], "mode": "bogus"})

    def test_species_rejects_non_string_mode(self):
        with pytest.raises(ValueError, match="mode must be a string"):
            TOOL_HANDLERS["species"]({"contigs": ["a.fna"], "mode": 5})

    def test_species_rejects_non_string_db_dir(self):
        with pytest.raises(ValueError, match="db_dir must be a string"):
            TOOL_HANDLERS["species"]({"contigs": ["a.fna"], "db_dir": 5})

    @pytest.mark.parametrize("bad", [0, -1, "4", 1.5, True])
    def test_species_rejects_bad_threads(self, bad):
        with pytest.raises(ValueError, match="threads must be an integer"):
            TOOL_HANDLERS["species"]({"contigs": ["a.fna"], "threads": bad})

    def test_db_setup_rejects_unknown_tier(self):
        with pytest.raises(ValueError, match="tier must be one of"):
            TOOL_HANDLERS["db_setup"]({"tier": "bogus"})

    def test_db_setup_rejects_missing_tier(self):
        with pytest.raises(ValueError, match="tier must be a non-empty string"):
            TOOL_HANDLERS["db_setup"]({})

    def test_db_setup_rejects_non_string_source(self):
        with pytest.raises(ValueError, match="source must be a string"):
            TOOL_HANDLERS["db_setup"]({"tier": "panel", "source": 7})

    def test_db_status_rejects_unknown_format(self):
        with pytest.raises(ValueError, match="format must be one of"):
            TOOL_HANDLERS["db_status"]({"format": "xml"})

    def test_db_list_rejects_unknown_format(self):
        with pytest.raises(ValueError, match="format must be one of"):
            TOOL_HANDLERS["db_list"]({"format": "xml"})


class TestEndToEnd:
    def test_species_call_over_the_wire(self, monkeypatch):
        monkeypatch.setattr(cli, "_identify_one", lambda c, m, d, threads=4: _payload(c))
        (response,) = exchange(_call("species", {"contigs": ["a.fna"], "mode": "marker"}))
        assert response["id"] == 1
        assert response["result"]["isError"] is False
        content = response["result"]["content"]
        assert content[0]["type"] == "text"
        assert json.loads(content[0]["text"])["verdict"]["species"] == "Salmonella"

    def test_handler_failure_renders_is_error(self, monkeypatch):
        def boom(contigs: str, mode: str, db_dir: str | None, threads: int = 4) -> dict[str, Any]:
            raise RuntimeError("kaboom")

        monkeypatch.setattr(cli, "_identify_one", boom)
        (response,) = exchange(_call("species", {"contigs": ["a.fna"]}))
        assert response["result"]["isError"] is True
        assert "RuntimeError: kaboom" in response["result"]["content"][0]["text"]

    def test_invalid_arguments_render_is_error(self):
        (response,) = exchange(_call("species", {"contigs": ["a.fna"], "mode": "bogus"}))
        assert response["result"]["isError"] is True
        assert "ValueError" in response["result"]["content"][0]["text"]


class TestSchemas:
    def test_roster_matches_handlers(self):
        assert [tool["name"] for tool in TOOLS] == list(TOOL_HANDLERS)

    def test_species_schema(self):
        (species,) = [tool for tool in TOOLS if tool["name"] == "species"]
        schema = species["inputSchema"]
        assert schema["required"] == ["contigs"]
        props = schema["properties"]
        assert props["contigs"] == {"type": "array", "items": {"type": "string"}}
        assert props["mode"]["enum"] == MODES
        assert props["mode"]["default"] == "marker"
        assert props["db_dir"] == {"type": "string"}
        assert props["threads"] == {"type": "integer", "minimum": 1, "default": 4}

    def test_db_status_schema(self):
        (tool,) = [tool for tool in TOOLS if tool["name"] == "db_status"]
        schema = tool["inputSchema"]
        assert schema["required"] == []
        assert schema["properties"]["format"]["enum"] == ["table", "tsv", "csv", "json"]
        assert schema["properties"]["format"]["default"] == "json"

    def test_db_list_schema(self):
        (tool,) = [tool for tool in TOOLS if tool["name"] == "db_list"]
        schema = tool["inputSchema"]
        assert schema["required"] == []
        assert schema["properties"]["format"]["enum"] == ["table", "tsv", "csv", "json"]
        assert schema["properties"]["format"]["default"] == "json"

    def test_db_setup_schema(self):
        (tool,) = [tool for tool in TOOLS if tool["name"] == "db_setup"]
        schema = tool["inputSchema"]
        assert schema["required"] == ["tier"]
        assert schema["properties"]["tier"]["enum"] == TIERS
        assert schema["properties"]["source"] == {"type": "string"}


class TestCliRegistration:
    def test_mcp_command_registered(self):
        assert "mcp" in main.commands

    def test_mcp_help_exits_zero(self):
        result = CliRunner().invoke(main, ["mcp", "--help"])
        assert result.exit_code == 0
        assert "MCP stdio server" in result.output

    def test_module_main_serves_stdin(self, monkeypatch, capsys):
        monkeypatch.setattr(
            "sys.stdin", StringIO('{"jsonrpc": "2.0", "id": 5, "method": "tools/list"}\n')
        )
        runpy.run_module("gside.mcp", run_name="__main__")
        response = json.loads(capsys.readouterr().out)
        assert response["id"] == 5
        assert response["result"]["tools"]


class TestMetaTool:
    def test_returns_meta_payload(self, monkeypatch):
        seen: dict[str, Any] = {}

        def fake_payload(reads, mode, db_dir, threads=4):
            seen.update(reads=reads, mode=mode, db_dir=db_dir, threads=threads)
            return {"analysis_type": "metagenomic_profiling", "verdict": {"richness": 0}}

        monkeypatch.setattr(cli, "_meta_payload", fake_payload)
        (response,) = exchange(_call("meta", {"reads": ["r1.fq", "r2.fq"], "mode": "sylph"}))
        payload = json.loads(response["result"]["content"][0]["text"])
        assert payload["analysis_type"] == "metagenomic_profiling"
        assert seen == {"reads": ["r1.fq", "r2.fq"], "mode": "sylph", "db_dir": None, "threads": 4}

    def test_rejects_unknown_mode(self):
        (response,) = exchange(_call("meta", {"reads": ["r.fq"], "mode": "bogus"}))
        assert response["result"]["isError"] is True
        assert "kraken2, sylph, all" in response["result"]["content"][0]["text"]

    def test_rejects_empty_reads(self):
        (response,) = exchange(_call("meta", {"reads": []}))
        assert response["result"]["isError"] is True

    def test_rejects_more_than_two_reads(self):
        (response,) = exchange(_call("meta", {"reads": ["a.fq", "b.fq", "c.fq"]}))
        assert response["result"]["isError"] is True


class TestMetaReadsValidation:
    def test_reads_must_be_array(self):
        (response,) = exchange(_call("meta", {"reads": "r.fq"}))
        assert response["result"]["isError"] is True
        assert "array of file path strings" in response["result"]["content"][0]["text"]
