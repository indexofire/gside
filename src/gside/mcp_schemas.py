"""MCP tools/list declarations: names, descriptions, inputSchemas (protocol: gside.mcp).

Pure data — the roster mirrors the handler map in :mod:`gside.mcp_tools`.
Descriptions are wire-visible wording; treat any edit like a schema change.
MODES/TIERS come from :mod:`gside.cli` so the wire enums can never drift from
the CLI surface.
"""

from __future__ import annotations

from gside.cli import META_MODES, MODES, TIERS

_FILES: dict[str, object] = {"type": "array", "items": {"type": "string"}}
_STR: dict[str, object] = {"type": "string"}
_FORMAT: dict[str, object] = {
    "type": "string",
    "enum": ["table", "tsv", "csv", "json"],
    "default": "json",
}


def _tool_entry(
    name: str,
    description: str,
    properties: dict[str, object] | None = None,
    required: list[str] | None = None,
) -> dict[str, object]:
    schema = {"type": "object", "properties": properties or {}, "required": required or []}
    return {"name": name, "description": description, "inputSchema": schema}


TOOLS: list[dict[str, object]] = [
    _tool_entry(
        "species",
        "Identify species from assembled contigs FASTA. Returns a JSON verdict with"
        " methods (marker/panel/mash_refseq/sourmash/gtdbtk), arbitration, and"
        " confidence.",
        dict(
            contigs=_FILES,
            mode={"type": "string", "enum": list(MODES), "default": "marker"},
            db_dir=_STR,
            threads={"type": "integer", "minimum": 1, "default": 4},
        ),
        ["contigs"],
    ),
    _tool_entry(
        "meta",
        "Profile a metagenomic sample (kraken2 and/or sylph) into a JSON species"
        " abundance verdict. Accepts FASTQ reads (single or paired) or contigs FASTA.",
        dict(
            reads=_FILES,
            mode={"type": "string", "enum": list(META_MODES), "default": "all"},
            db_dir=_STR,
            threads={"type": "integer", "minimum": 1, "default": 4},
        ),
        ["reads"],
    ),
    _tool_entry(
        "db_status",
        "Check database readiness status. Returns JSON with per-tier ready/path/tier info.",
        dict(format=_FORMAT),
    ),
    _tool_entry(
        "db_list",
        "List available database tiers with readiness. Returns JSON with tier metadata.",
        dict(format=_FORMAT),
    ),
    _tool_entry(
        "db_setup",
        "Install or update a database tier. Progress on stderr; result JSON on return.",
        dict(tier={"type": "string", "enum": list(TIERS)}, source=_STR),
        ["tier"],
    ),
]
