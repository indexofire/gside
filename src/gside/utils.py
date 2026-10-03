"""Shared parsing helpers (ported verbatim from hermes-bacmap utils)."""

from __future__ import annotations


def parse_db_header(sseqid: str) -> tuple[str, str, str, str]:
    """Parse a 'db~~~gene~~~accession~~~product' FASTA/BLAST header.

    Returns (gene, accession, product, reserved); missing fields are "".
    """
    fields = sseqid.split("~~~")
    if len(fields) >= 4:
        return fields[1].strip(), fields[2].strip(), fields[3].strip(), ""
    if len(fields) >= 3:
        return fields[1].strip(), fields[2].strip(), "", ""
    if len(fields) >= 2:
        return fields[1].strip(), "", "", ""
    return sseqid.strip(), "", "", ""
