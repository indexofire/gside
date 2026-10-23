"""Canonical species labels for the JSON contract.

bacmap (the downstream consumer) matches species labels with
``[A-Za-z]+(?:_[A-Za-z0-9]+)+`` and derives three-letter species codes from the
first two tokens, so every tier must emit ``Genus_species``:

* panel / mash_refseq — ``Streptococcus agalactiae 2603V/R`` → ``Streptococcus_agalactiae``
* sourmash — ``...;s__Vibrio parahaemolyticus`` → ``Vibrio_parahaemolyticus``
* gtdbtk — ``Vibrio cholerae`` → ``Vibrio_cholerae``
* marker — already registry labels (``Salmonella_enterica``, ``Vibrio_parahaemolyticus``,
  ``Escherichia_coli``, ``Shigella``), emitted as-is.

Strain / accession detail lives in a separate ``strain`` field; it is never
glued into ``species`` (a dot in the label breaks bacmap's matcher).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

UNKNOWN = "Unknown"

_RANK_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")
_TOKEN_RE = re.compile(r"[^._\s/]+")
# sourmash names look like "GCF_001457695.1 Campylobacter jejuni".
_ACCESSION_PREFIX_RE = re.compile(r"^GC[AF]_\d+(?:\.\d+)?\s+")
_BARE_ACCESSION_RE = re.compile(r"GC[AF]_\d+(?:\.\d+)?")

# Empty / placeholder verdicts belong to no species.
_PLACEHOLDERS = frozenset(
    {
        "",
        "unknown",
        "mixed",
        "mixed/unknown",
        "n/a",
        "na",
        "none",
        "-",
        "unidentified",
        "unclassified",
        "pending_gside",
    }
)

# `Campylobacter sp.`, `Listeria spp.`, `Uncultured bacterium` are not species.
_NON_SPECIES_TOKENS = frozenset({"sp", "spp", "species", "bacterium", "uncultured", "cf"})
_NON_SPECIES_RE = re.compile(r"sp\d+")
# GTDB disambiguation suffixes: "Campylobacter_D jejuni", "Vibrio_A cholerae".
_GTDB_SUFFIX_RE = re.compile(r"[A-Z]")


@dataclass(frozen=True)
class SpeciesName:
    label: str
    display: str
    strain: str


_UNKNOWN_NAME = SpeciesName(label=UNKNOWN, display=UNKNOWN, strain="")


def _strip_rank(text: str) -> str:
    for prefix in _RANK_PREFIXES:
        if text.startswith(prefix):
            return text[len(prefix) :].strip()
    return text


def _species_from_lineage(text: str) -> str | None:
    """Last ``s__`` rank of a GTDB lineage; None when the lineage has none."""
    for part in reversed(text.split(";")):
        part = part.strip()
        if part.startswith("s__"):
            return part[3:].strip()
    return None


def parse_species_name(raw: str) -> SpeciesName:
    """Split an arbitrary tier label into (canonical label, display, strain)."""
    text = (raw or "").strip()
    if ";" in text:
        lineage_species = _species_from_lineage(text)
        if lineage_species is None:
            return _UNKNOWN_NAME
        text = lineage_species
    text = _strip_rank(_ACCESSION_PREFIX_RE.sub("", text))
    if text.lower() in _PLACEHOLDERS or _BARE_ACCESSION_RE.fullmatch(text):
        return _UNKNOWN_NAME

    matches = list(_TOKEN_RE.finditer(text))
    if not matches:
        return _UNKNOWN_NAME
    tokens = [m.group(0) for m in matches]
    keep = list(range(len(tokens)))
    # Drop a GTDB single-letter suffix on the genus or the species epithet.
    if len(keep) >= 3 and _GTDB_SUFFIX_RE.fullmatch(tokens[keep[1]]):
        keep.pop(1)
    if len(keep) >= 3 and _GTDB_SUFFIX_RE.fullmatch(tokens[keep[2]]):
        keep.pop(2)
    kept = [tokens[i] for i in keep]

    second = kept[1].lower() if len(kept) >= 2 else ""
    if second in _NON_SPECIES_TOKENS or _NON_SPECIES_RE.fullmatch(second):
        return _UNKNOWN_NAME

    strain = text[matches[keep[2]].start() :].strip(" ._/-") if len(keep) > 2 else ""
    return SpeciesName(
        label="_".join(kept[:2]),
        display=" ".join(kept[:2]),
        strain=strain,
    )


def canonical_species(raw: str) -> str:
    """``raw`` → ``Genus_species`` (``Unknown`` for placeholders)."""
    return parse_species_name(raw).label


def display_name(label: str) -> str:
    """Expand an already-canonical label for humans (``Vibrio_cholerae_O1`` →
    ``Vibrio cholerae O1``). Never reduces, so serotype suffixes survive."""
    text = (label or "").strip()
    return text.replace("_", " ") if text else UNKNOWN


def apply_species_name(result: dict[str, Any], raw: str) -> dict[str, Any]:
    """Write canonical ``species`` (+ ``species_raw`` / ``strain``) into ``result``."""
    parsed = parse_species_name(raw)
    result["species"] = parsed.label
    if raw.strip() and raw.strip() != parsed.label:
        result["species_raw"] = raw
    if parsed.strain:
        result["strain"] = parsed.strain
    return result
