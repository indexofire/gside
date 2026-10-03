# AGENTS.md

Agent guidance for this repository. Read this before writing or editing any code.

---

## Project Overview

- **Language**: Python 3.12 (standalone species-identification CLI, no frontend)
- **Package manager**: [pixi](https://pixi.sh) — conda tools (blast/skani/mash/mmseqs2) + PyPI editable install via `pixi.toml`
- **Linter + Formatter**: [ruff](https://docs.astral.sh/ruff/) (config in `pyproject.toml`; run via pip-installed dev extras, no pixi task)
- **Testing**: pytest (`tests/`; single file `tests/unit/test_cli.py`, needs real `blastn`)
- **Docs**: MkDocs Material + static-i18n, English default + Chinese secondary (`docs/*.md` + `*.zh.md` pairs)

---

## Environment Setup

```bash
# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Install conda tools + editable package
pixi install
pip install -e .

# Optional docs extras (dev environment)
pip install -e .[dev]
```

Bioinformatics binaries come from the pixi `default` env; the Python package is
pip-installed (editable). `pip install -e .[dev]` adds pytest/ruff/mypy/mkdocs.

---

## Build / Run Commands

```bash
# Species identification (5 modes, default marker)
pixi run -e dev python -m gside species contigs.fna --mode marker
pixi run -e dev python -m gside species contigs.fna --mode all

# Database management
pixi run -e dev python -m gside db status
pixi run -e dev python -m gside db setup --tier panel   # panel (default) | mash | mini | all

# After pip install -e ., the console script is also available
gside species contigs.fna --mode panel
gside --version
```

Check `pixi.toml` for the canonical task definitions — only `docs` / `docs-build`
exist as pixi tasks; everything else runs via `python -m` or the `gside` script.

---

## Lint / Format Commands

```bash
# No pixi tasks exist for lint; use the pip-installed ruff from .[dev]
python -m ruff check src/ tests/
python -m ruff format --check src/ tests/
python -m ruff format src/ tests/
```

Config: `pyproject.toml` (`target-version = "py312"`, `line-length = 100`).
**Always run ruff before committing.** Fix all lint errors; do not use `# noqa`
suppression unless unavoidable, and always add a comment explaining why.

---

## Testing Commands

```bash
# Run all tests (dev env now ships pytest + pytest-cov)
pixi run -e dev pytest tests/ -q

# Run with verbose output
pixi run -e dev pytest tests/ -v

# Unit tests only (pure logic, no binaries/databases — always fast)
pixi run -e dev pytest tests/unit -q

# Integration tests (need blastn/skani/mash + databases)
pixi run -e dev pytest tests/integration -q

# With coverage (ratchet: must not drop below the committed floor)
pixi run -e dev pytest tests/ -q --cov=src/gside --cov-report=term-missing
```

Layout: `tests/unit/` = pure, hermetic, milliseconds; `tests/integration/` =
needs real binaries or databases (gated by `skipif`, never mock a binary —
either run it or skip). `TestSpeciesMarker` needs `blastn` resolvable via
`GSIDE_PIXI_BIN` or pixi env (discovery order: `src/gside/config.py`).
CI (`.github/workflows/test.yml`, pip channel): ruff check + format --check,
`pyright src/gside/` (must stay zero), `pytest -q --cov=src/gside`.
Binary-needing tests skip gracefully without conda tools; full runs happen
in the pixi dev env locally.
Baseline (2026-10-03): 330 tests, total coverage **97%** — all modules 100%
except: ani 93%, multigene 98%, taxonomic 93% (quarantined), cli 99%,
db 92%, kma 94%, kmer 99%, read_mapper 99%. Enforced by `[tool.coverage.report] fail_under = 95` in `pyproject.toml`;
raise it every time coverage climbs.

---

## Test-Driven Development (TDD)

TDD is mandatory for all behavior changes and bug fixes. No exceptions.

### The loop

1. **Red** — write a failing test FIRST that pins the desired behavior
   (unit test if the logic is pure, integration test if it needs a binary/DB).
   Run it, watch it fail for the right reason.
2. **Green** — write the minimal implementation that makes it pass.
3. **Refactor** — clean up with tests green; re-run the full suite + ruff.

### Rules

- **One behavior per test** — a test name states the behavior
  (`test_ani_beats_marker`, not `test_arbitrate2`).
- **New/changed code ships with tests.** A PR without a covering test is
  incomplete: pure helpers → `tests/unit/`; binary/DB paths →
  `tests/integration/` with `skipif` guards, never mocked binaries.
- **Bug fix = regression test first.** Reproduce with a failing test, then fix.
- **Coverage ratchet.** Canonical command is `pixi run -e dev cov` (explicit
  `--cov-fail-under`). The `fail_under` value in `pyproject.toml` IS enforced
  by pytest-cov whenever `--cov` runs without an explicit flag value, so the
  config and the flags must stay equal. CI checkouts lack the gitignored
  databases, so the CI ceiling sits ~1pt below local (CI ~96% vs local 97%)
  — the floor is 95 to hold in BOTH environments; local coverage must still
  never go below its last high-water mark. Raise the floor in `pyproject.toml`,
  the `cov` task, the CI workflow, and this file together — `test_floors_agree`
  fails otherwise. Capture exit codes **without pipes** (`cmd > log; echo $?`),
  pipes mask them.
- **Deterministic tests only.** No network, no wall-clock, no absolute paths
  outside the repo (use `tmp_path` + `monkeypatch`). Real binaries are fine;
  the network is not.
- **Fast unit suite.** `tests/unit` must finish in seconds; anything slower
  belongs in `tests/integration`.

---

## Code Style Guidelines

### Formatting
- Line length: **100 characters** (see `pyproject.toml`; not the ruff default 88)
- Indentation: **4 spaces** (never tabs)
- String quotes: **double quotes** (`"..."`)
- Keep functions small; JSON payloads are built as plain dicts

### Imports
- Use **absolute imports** only (`from gside.config import ...`); no relative imports
- Import order: stdlib → third-party → local (ruff isort-compatible)
- No wildcard imports (`from module import *`)

### Naming Conventions
| Construct | Convention | Example |
|---|---|---|
| Module/file | `snake_case` | `multigene_identifier.py` |
| Function | `snake_case` | `identify_multigene()` |
| Variable | `snake_case` | `db_dir` |
| Class | `PascalCase` | `MultiGeneResult` |
| Constant | `UPPER_SNAKE_CASE` | `_MIN_IDENTITY` |
| Private | leading underscore | `_arbitrate()` |

### Type Annotations
- **Always annotate** function signatures (parameters + return type)
- Use built-in generics (`list[str]`, `dict[str, Any]`) over `typing.List` / `typing.Dict`
- Use `X | Y` union syntax over `Union[X, Y]`

### Error Handling
- `gside species` **always exits 0** — per-method failures go into
  `methods.<mode>.error`, never raise through the CLI
- `gside db setup` returns 1 only when a result message contains `"ERROR"`
- Never silence exceptions with empty `except` blocks
- Prefer `pathlib.Path` over `os.path` for filesystem operations

### CLI Patterns
- Click group style (mirrors gmlst): `HELP_SETTINGS` (`-h`/`--help`), `version_option`, global `-v`/`-q`, thin commands registered via `add_command`
- Subcommands: `species` / `db` group (`status|setup|list`) / `validate` / `--version`; entry `gside.cli:main`
- `--mode` choices are exactly `marker|panel|mash_refseq|sourmash|all` — do not invent modes
- `--db-dir` overrides the database root for one invocation

### File & Module Organization
```
src/gside/
├── __init__.py
├── cli.py                    # Entry point; all-mode runner + _arbitrate()
├── config.py                 # DATA_DIR / SPECIES_DB_DIR / pixi binary discovery
├── db.py                     # gside db status/setup/list + downloaders
├── engine/                   # Algorithm abstraction (mostly unused, see gotchas)
│   ├── hits.py               #   Unified Hit dataclass
│   ├── registry.py           #   Backend registry
│   ├── read_mapper.py        #   ReadMapper facade
│   └── backends/             #   blast / skani / kmer / minimap2 / kma / mmseqs2
└── analysis/
    ├── multigene_identifier.py   # L1 marker rules (blastn, live)
    ├── ani_identifier.py         # L2 panel/mash (live)
    ├── sourmash_identifier.py    # L2 sourmash (needs external binary + GTDB db)
    ├── taxonomic_validator.py    # L3 (gside validate; standard needs CheckM2/GTDB-Tk env)
    └── species_canon.py          # DORMANT (no importers)
tests/
└── unit/
    └── test_cli.py           # Real-BLAST smoke + arbitration unit tests
data/
├── db/
│   ├── D1_marker/            # mini tier, TRACKED in git: rules + fasta + blastdb
│   ├── D2_ani/               # panel tier (skani sketch + genomes), gitignored
│   │                         # except manifests/ (tracked build recipe)
│   ├── D3_mash/              # mash tier (mash.msh), gitignored
│   └── D4_sourmash/          # GTDB reps k31 + lineages, via db setup (3.9GB)
docs/
├── index.md / index.zh.md    # + architecture|installation|usage|reference/, same pairing
└── (nav in mkdocs.yml; en default, zh secondary; build with --strict)
```

- Keep modules focused and small; no circular imports
- L1 runtime data lives in `data/db/D1_marker/` (tracked); large L2 databases are
  gitignored and provisioned via `gside db setup`

---

## pyproject.toml Conventions

```toml
[tool.ruff]
target-version = "py312"
line-length = 100

[tool.pytest.ini_options]
testpaths = ["tests"]
```

---

## Key Rules for Agents

1. **Run `pixi run -e dev check` after every code change** (ruff check + format check; `fix` to auto-fix)
2. **Never use `# type: ignore` or `# noqa` without an explanatory comment**
3. **Never commit secrets** — tunables go into `marker_rules.yaml` / env vars, not code
4. **Prefer `pathlib.Path` over `os.path`** for all filesystem operations
5. **All functions must have type annotations** — no untyped signatures
6. **One responsibility per module** — don't let files grow into god-modules
7. **Don't add conda deps without updating `pixi.toml`** — use `pixi add <pkg>`; PyPI-only packages belong in `[feature.dev.pypi-dependencies]`, never in `[dependencies]` (conda channels don't carry them, e.g. `mkdocs-static-i18n`)
8. **Do not rename `species_markers_v2`** — it is the JSON output contract name (`database.name`), not a filename
9. **Do not wire up `species_canon.py`** — dormant with no importers; touching it is a separate project (`taxonomic_validator.py` is wired via `gside validate` since Cycle 19)
10. **Docs come in en/zh pairs** — every `docs/**/*.md` change needs the matching `.zh.md` update, and `mkdocs build --strict` must pass
11. **TDD is non-negotiable** — failing test first, then implementation; every new line covered; coverage never goes down (see TDD section)
12. **Keep `pyright` at zero and CI green** — `pixi run -e dev typecheck` must report 0 errors (quarantined files are excluded in `pyproject.toml`, not waived); `.github/workflows/test.yml` must pass before merge

## Known Gotchas (read before touching these areas)

- **Binary discovery**: `GSIDE_PIXI_BIN` env → repo `.pixi` → `~/.pixi` → hardcoded hermes-bacmap fallback path → `PATH` (`config.py`). The hardcoded path is machine-specific; prefer `GSIDE_PIXI_BIN` on other machines.
- **`gside species` never fails loud**: check `methods.<mode>.error` in JSON, not the exit code.
- **Marker thresholds are global**: `_MIN_IDENTITY=85.0` / `_MIN_COVERAGE=60.0` / `_HIGH_CONF=90.0`; rule-level `min_coverage` in YAML is currently **ignored** by the engine; `priority_over` is dead code.
- **Panel thresholds**: ANI ≥ 95 + AF ≥ 0.65 → high; 93–95 → medium. Mash identity ≥ 0.97 → high, 0.90–0.97 → medium.
- **Arbitration is layer-first**: `(layer, confidence)` tuple — any L2 hit beats any L1 hit (`cli.py:_arbitrate`).
- **`sourmash` needs its GTDB db**: provisioned via pixi (`sourmash>=4.9`) + `gside db setup --tier sourmash` (~3.9GB rs226 reps, renamed in place); `gather` output is read from `-o` file, never `--csv -` (unsupported flag).
- **DB download verification**: `db.py` pins digests (`MASH_MD5`, `PANEL_SHA256`, `SOURMASH_LINEAGES_SHA256`); empty constants (`SOURMASH_SIG_SHA256`) mean warn-and-skip. Tar extraction uses `filter="data"`. All pinned artifacts live in the GitHub release tag `db-v0.1` (panel) — re-pin `PANEL_SHA256` whenever the release asset is replaced.
- **`engine/` is wired for the species-ID hot paths**: marker → `BlastBackend`, panel → `SkaniBackend`, mash_refseq → `MashBackend` (analysis layer resolves binaries via `require_bin` and injects them with `binary=`; outputs are byte-equivalent to the former inline subprocess calls). Shared scaffolding lives in `engine/_env.py` (`require_bin`/`run_checked`; backends pass their module-level `which` as `resolver=` to keep test monkeypatch points). Still unconsumed: `KmaBackend`, `Mmseqs2Backend`, `MinimapBackend`, `ReadMapper`, `SequenceMatcher` — `KmaBackend` imports cleanly (`parse_db_header` ported verbatim into `gside/utils.py`) but has no in-repo callers and no provisioned binary.
- **FASTA header contract**: `>markers~~~<gene>~~~<desc> role=<primary|confirm|typing|virulence>`; parser takes `split("~~~")[1]` — the `markers` prefix itself is never validated.
