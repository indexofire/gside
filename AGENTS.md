# AGENTS.md

Agent guidance for this repository. Read this before writing or editing any code.

---

## Project Overview

- **Language**: Python 3.12 (standalone species-identification CLI, no frontend)
- **Package manager**: [pixi](https://pixi.sh) — conda tools (blast/skani/mash/sourmash/kraken2/bracken/sylph) + PyPI editable install via `pixi.toml`
- **Linter + Formatter**: [ruff](https://docs.astral.sh/ruff/) (config in `pyproject.toml`; `pixi run -e dev check`, or the pip-installed dev extras)
- **Testing**: pytest — `tests/unit/` hermetic (no binaries), `tests/integration/` needs real binaries/databases (skipif-guarded)
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
Baseline (2026-10-21): **572 tests**, total coverage **98.9%**
— all modules 100% except: ani 93%, multigene 99%, sourmash 99%,
taxonomic 93% (quarantined), db 99%, kraken2 99%, progress 97%.
Enforced by `[tool.coverage.report] fail_under = 95` in `pyproject.toml`;
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
- Click group style (mirrors gmlst/gapit family): `HELP_SETTINGS` (`-h`/`--help`), bare invocation prints help and exits 2 (`no_args_is_help`), `-v`/`-q` mutually exclusive (UsageError), thin commands registered via `add_command`
- Every option carries a short flag: species `-m/--mode -f/--format -D/--db-dir -o/--output -j/--jobs -t/--threads -n/--no-header -a/--auto-setup`; meta `-m/--mode -f/--format -D/--db-dir -o/--output -s/--samples -j/--jobs -t/--threads -n/--no-header -a/--auto-setup`; db `-t/--tier -s/--source -f/--format`
- `-D/--db-dir`, `-f/--format` and `-t/--threads` are built by the shared factories `_db_dir_option` / `_format_option` / `_threads_option` in `cli.py` — add new commands through them so flags never drift
- Agent-friendly outputs: species defaults to JSON on stdout; `-o` writes the payload to a file with gapit-style streaming (tsv/md rows flushed per input as they complete — a mid-batch failure persists the prefix; json single document at the end; stdout stays empty, notice on stderr); `db status|list` accept `--format table|tsv|json`; `db setup` accepts `--format text|json` (json routes progress to stderr); `-V -J` prints version as JSON
- Subcommands: `species` / `meta` / `db` group (`status|setup|list`) / `validate` / `mcp` (MCP stdio server) / `--version`; entry `gside.cli:main`
- `meta` batch mode: `-s/--samples sample_id<TAB>r1[<TAB>r2]` + `-j/--jobs` (order preserved); JSON becomes an array with a `sample` key, tsv/csv gain a leading `sample` column. `--samples` and positional reads are mutually exclusive (UsageError)
- `species --mode` choices are exactly `marker|panel|mash_refseq|sourmash|gtdbtk|all`; `meta --mode` choices are exactly `kraken2|sylph|all` — do not invent modes
- `--db-dir` overrides the database root for one invocation
- `verdict` shape is shared by single-mode and `all`: `species` (canonical), `display_name`, `strain`, `confidence`, `basis` (winner only), `supporting`, `agreement`

### File & Module Organization
```
src/gside/
├── __init__.py
├── cli.py                    # Entry point; all-mode runner + _arbitrate()
├── config.py                 # DATA_DIR / SPECIES_DB_DIR / pixi binary discovery
├── db.py                     # gside db status/setup/list + downloaders
├── mcp.py                    # MCP stdio server: JSON-RPC protocol + serve loop (`gside mcp`)
├── mcp_schemas.py            # MCP tools/list declarations (species/meta/db_status/db_list/db_setup)
├── mcp_tools.py              # MCP tool handlers (call cli._identify_one + db callables)
├── engine/                   # Alignment backends actually used by analysis/
│   ├── _env.py               #   require_bin / run_checked / pixi_path
│   ├── hits.py               #   BLAST tabular Hit row
│   └── backends/             #   blast / skani / kmer (mash + sourmash)
└── analysis/
    ├── _envelope.py              # Shared {analysis_type,method,database,result}
    ├── _thresholds.py            # Shared ANI_HIGH/ANI_MEDIUM + ani_confidence()
    ├── species_name.py           # Genus_species canonicalizer (bacmap contract)
    ├── multigene_identifier.py   # L1 marker rules (blastn, live)
    ├── ani_identifier.py         # L2 panel/mash (live)
    ├── sourmash_identifier.py    # L2 sourmash (needs external binary + GTDB db)
    ├── gtdbtk_identifier.py      # L3 GTDB-Tk placement (needs external gtdbtk + GTDBTK_DATA_PATH)
    ├── kraken2_identifier.py     # meta: kraken2 + bracken
    ├── sylph_identifier.py       # meta: sylph ANI profiling (FASTQ only)
    └── taxonomic_validator.py    # L3 (gside validate; standard needs CheckM2/GTDB-Tk env)
tests/
├── unit/                     # Hermetic, no binaries, seconds
│   └── test_cli.py           # Rendering + arbitration only (pure)
└── integration/              # Real binaries/databases, skipif-guarded
    └── test_marker_smoke.py  # Real blastn marker calls
data/
├── db/
│   ├── D1_marker/            # mini tier, TRACKED in git: rules + fasta + blastdb
│   ├── D2_ani/               # panel tier (skani sketch + genomes), gitignored
│   │                         # except manifests/ (tracked build recipe)
│   ├── D3_mash/              # mash tier (mash.msh), gitignored
│   ├── D4_sourmash/          # GTDB reps k31 + lineages, via db setup (3.9GB)
│   ├── D5_kraken2/           # meta tier (Kraken2 Standard-8), manual, gitignored
│   └── D6_sylph/             # meta_sylph tier (.syldb), manual, gitignored
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
9. **Keep `engine/` limited to backends the analysis layer calls** — blast / skani / kmer only. The old name-registry, ReadMapper, minimap2, kma and mmseqs2 abstractions were deleted as unconsumed ports from gmlst; do not reintroduce scaffolding without an in-repo caller and a test.
10. **Docs come in en/zh pairs** — every `docs/**/*.md` change needs the matching `.zh.md` update, and `mkdocs build --strict` must pass
11. **TDD is non-negotiable** — failing test first, then implementation; every new line covered; coverage never goes down (see TDD section)
12. **Keep `pyright` at zero and CI green** — `pixi run -e dev typecheck` must report 0 errors (quarantined files are excluded in `pyproject.toml`, not waived); `.github/workflows/test.yml` must pass before merge

## Known Gotchas (read before touching these areas)

- **Binary discovery**: `GSIDE_PIXI_BIN` env → repo `.pixi/envs/default/bin` → `~/.pixi/envs/default/bin` → `PATH` (`config.py`). There is no machine-specific path in the package; point `GSIDE_PIXI_BIN` at a conda env when the tools live elsewhere.
- **`gside species` never fails loud**: check `methods.<mode>.error` in JSON, not the exit code.
- **Marker thresholds**: global `_MIN_IDENTITY=85.0` / `_MIN_COVERAGE=60.0` / `_HIGH_CONF=90.0` gate `detected_markers` and `exclude_genes`. Rule matching uses each rule's own `min_identity`/`min_coverage` (defaults 90/60) against every raw hit of the gene, so the two *Campylobacter* rules with `min_coverage: 50` really do accept 50-60% coverage. `priority_over` is supported (priority species suppress non-priority ones) but no shipped rule uses it.
- **Panel thresholds**: ANI ≥ 95 + AF ≥ 0.65 → high; 93–95 → medium. Mash identity ≥ 0.97 → high, 0.90–0.97 → medium.
- **Arbitration is layer-first**: `(layer, confidence-rank)` tuple — gtdbtk=3 beats any L2 hit, any L2 hit beats any L1 hit, and within a layer `medium` now beats `low` (`cli.py:_arbitrate`). `verdict.basis` stays the single winning method (bacmap maps it to a validated `method`); `verdict.supporting` + `verdict.agreement` expose cross-method corroboration (exact label match — `DEC` vs `Escherichia_coli` is intentionally a non-match; do not alias them).
- **Missing-database errors are structured**: alongside `methods.<mode>.error`, the payload carries `suggested_command` / `missing_tier` / `estimated_size` and `help_url` (from `gside.db._TIERS`). Both `species --auto-setup` and `meta --auto-setup` install the tier once (progress on stderr) and retry that method; meta's batch mode installs on the first sample, then fans out. With an explicit `--db-dir`, auto-setup is skipped (user-managed root). GTDB-Tk's own database is external (`GTDBTK_DATA_PATH`/`GTDBDB`) — it is NOT a `db setup` tier.
- **`sourmash` needs its GTDB db**: provisioned via pixi (`sourmash>=4.9`) + `gside db setup --tier sourmash` (~3.9GB rs226 reps, renamed in place); `gather` output is read from `-o` file, never `--csv -` (unsupported flag). Confidence comes from the gather containment ANI (`average_containment_ani`, 95/93 bands, shared with panel); the old `f_unique_weighted` 0.90/0.70 bands remain only as a fallback when that column is absent. `possible_mixture` requires >=2 **distinct species** at >=10% with ANI >=93 — multiple strains of one species do not trigger it.
- **Metagenomic tiers auto-download**: `db setup --tier meta` pulls Ben Langmead's Standard-8 (5.5GB, includes Bracken 50-300mer) and verifies each extracted file against the published md5 list; readiness needs `database150mers.kmer_distrib`. `--tier meta_sylph` pulls the 14GB two-stage GTDB r226 `.syl2db` (MD5 pinned). sylph accepts both `.syldb` and `.syl2db`; `.syl2db` is preferred and needs sylph >= 1.0.
- **DB download verification**: panel pins SHA256, mash MD5 (Zenodo only publishes MD5), sourmash lineages + signature SHA256, sylph DB MD5 (GCS `x-goog-hash`). Kraken2 is verified per extracted file against the published `.md5` manifest. The sourmash signature zip has no upstream digest; `SOURMASH_SIG_SHA256` was pinned after a locally verified download (complete size, zip integrity, 143,385 k=31 DNA signatures, CC0) and `GSIDE_SOURMASH_SIG_SHA256` overrides it at setup time. Tar extraction uses `filter="data"`. All pinned artifacts live in the GitHub release tag `db-v0.1` (panel) — re-pin `PANEL_SHA256` whenever the release asset is replaced.
- **`engine/` holds only live backends**: marker → `BlastBackend`, panel → `SkaniBackend`, mash_refseq/sourmash → `MashBackend`/`SourmashBackend` (analysis resolves binaries via `require_bin` and injects them with `binary=`; backends pass their module-level `which` as `resolver=` so tests can monkeypatch discovery). Shared scaffolding is `engine/_env.py` (`require_bin`/`run_checked`).
- **Species name contract is `Genus_species`**: every tier emits the underscored label via `analysis/species_name.py` (`apply_species_name`) and keeps the tool string in `species_raw`, the strain in `strain`. The normalizer strips GTDB `s__` ranks, a leading `GCF_/GCA_` accession (sourmash names look like `GCF_... Campylobacter jejuni`), and GTDB disambiguation suffixes (`Campylobacter_D jejuni` → `Campylobacter_jejuni`, `Vibrio_A cholerae` → `Vibrio_cholerae`); `sp.`/`sp036004205`/`unidentified`/`Unknown` collapse to `Unknown`. marker labels are already contract form (`Salmonella`, `V_parahaemolyticus`, `Vibrio_cholerae_O1`, `Shigella_EIEC`, `DEC`) and pass through untouched — never "reduce" them to two tokens. bacmap's `canonical_from_verdict` only accepts `[A-Za-z]+(?:_[A-Za-z0-9]+)+`, so a spaced or dotted label breaks it.
- **FASTA header contract**: `>markers~~~<gene>~~~<desc> role=<primary|confirm|typing|virulence>`; parser takes `split("~~~")[1]` — the `markers` prefix itself is never validated.
