# gseek engine migration

`gside` can route its marker BLAST/FASTQ screening, mash distance, sourmash
gather, and panel/skani lookup through
[`gseek`](https://github.com/indexofire/gseek) behind an explicit switch. The
default remains the existing subprocess implementation; gseek is optional.

## Switch

```bash
GSIDE_ENGINE=legacy   # default: existing backends, no gseek import
GSIDE_ENGINE=gseek    # route all migrated paths through gseek
GSIDE_ENGINE=shadow   # run both, log differences, return the legacy result
```

`GSIDE_ENGINE` is read per call, so a process can be tested with either mode
without restarting. Invalid values fail before any tool runs. The selector,
shadow comparison, and fallback are provided by `gseek.compat.Engine`; legacy
mode still avoids importing gseek entirely.

## Legacy is retained

This is a deliberate dual-path design, not a temporary migration leftover:

- `legacy` is the default, so installations without gseek behave exactly as
  before;
- every migrated path keeps its legacy implementation and remains usable;
- `shadow` runs both implementations, logs coarse differences, and returns the
  legacy result;
- legacy is the reference for parity tests and the fallback if gseek is
  unavailable;
- modes without a gseek adapter (`gtdbtk`, `kraken2`, `sylph`) are legacy-only.

Install gseek for local testing:

```bash
pixi install -e dev            # dev feature includes gseek>=0.10.0
# or: pip install -e ".[dev]"
# or, without installing a release:
PYTHONPATH=../gseek/src pixi run -e dev pytest tests/integration/test_gseek_parity.py -q
```

> The mash/sourmash paths require gseek 0.10.0 (`external_reference`), which is
> published on PyPI: `pip install gseek==0.10.0`. The optional extra is pinned
> accordingly. Development/test environments install it via the `dev` extra
> (`pip install -e ".[dev]"`) and the pixi `dev` feature; the sibling checkout
> remains useful for working on gseek itself.

## Migrated paths

| Path | Legacy | gseek |
| --- | --- | --- |
| Marker BLAST | `BlastBackend` + `markers_blastdb` | `gseek.align(backend="blast", target=markers.fasta)` |
| FASTQ marker screening | minimap2 PAF + local aggregation | `gseek.align(backend="minimap2", ...)` + same aggregation |
| mash_refseq (assembly/reads) | `MashBackend` against `mash.msh` | `gseek.distance(backend="mash", external_reference=True)` |
| sourmash gather | sketch/gather CSV parsing | `gseek.distance(backend="sourmash", external_reference=True)`, one run per query file |
| panel (skani) | `SkaniBackend` + `panel.sketch` DB | `gseek.distance(backend="skani", external_reference=True)` against the pre-built DB directory |

The gseek path deliberately preserves gside semantics:

- coverage is still `aligned_length / length` (not gseek span coverage);
- screening identity is still `nmatch / alen`, via `record.matches`;
- BLAST uses `task="megablast"`, `word_size=11`, and dust on to match the
  historical `blastn` defaults;
- tool discovery uses the same `GSIDE_PIXI_BIN`/pixi/PATH order.

## Not migrated yet

All gseek-backed gside modes are migrated (marker, mash, sourmash, panel).
Modes without a gseek adapter (`gtdbtk`, `kraken2`, `sylph`) remain on their
existing implementations.

## Verification

- `tests/unit/test_gseek_bridge_units.py`: record mapping, aggregation, mode
  dispatch, and fake-gseek parameter wiring.
- `tests/integration/test_gseek_parity.py`: real marker DB/genome comparison
  against the legacy path, plus FASTQ screening parity. Skipped when gseek is
  not importable.

The same gates apply on `main` (`ruff`, `pyright`, `pytest tests/unit --cov ...`,
plus integration). The parity test is collected by the default test suite when
gseek is installed; cases that need external binaries skip when those binaries
are absent.


## Local verification (2026-10-24)

`gside scan -m marker -f json` was run against the tracked N16961 genome and a
synthetic FASTQ generated from the `inva` marker:

| Input | legacy | gseek | shadow |
| --- | --- | --- | --- |
| `GCF_000006745.1` assembly | `Vibrio_cholerae` / high / `ompw` 94.767 / 100.1 | identical | identical |
| 21 reads covering `inva` | `Salmonella_enterica` / high / `inva` 100.0 / 98.3 | identical | identical |
| K12 assembly, `-m mash_refseq` | `Escherichia_coli` / high / identity 1.0 | identical | identical |

`shadow` produced no difference warnings. Unit suite: 748 passed; integration: 28 passed, 1 skipped
(with gseek installed; 6 are parity, 22 without gseek); coverage 96.88%.



## Real-data spot check (`~/temps/test`, 2026-10-24)

Inputs: one RefSeq genome from each of *Acinetobacter baumannii*,
*Campylobacter coli*, *Escherichia coli*, *Salmonella enterica*, and
*Shigella flexneri*. Each was run with
`gside scan <genome> -m marker+mash_refseq -f json` under all three engines.

- legacy, gseek, and shadow JSON were **byte-identical** for all five genomes.
- Example calls remained identical: `A. baumannii` marker `rpob-aba` 99.068 /
  `oxa51` 97.39, mash 0.978914; `S. flexneri` marker `ipah` 100.0, mash
  0.993922.
- A 260,000-read (150 bp, ~10x) FASTQ from *A. baumannii* produced identical
  marker (`Acinetobacter_baumannii`, high, `rpob-aba` 94.84 / `oxa51` 93.56)
  and mash (`0.9789136`, high) results in all three modes.
- panel (skani) on *A. baumannii* and *C. coli*: legacy/gseek/shadow JSON were
  byte-identical (A. baumannii 98.06 ANI / 85.79 AF; C. coli 99.72 / 98.09).
- A mini sourmash signature DB built from three RefSeq genomes produced
  identical `species`, `confidence`, `f_unique_weighted`, `ani`, and
  `gather_partition` for *A. baumannii* and *S. enterica* in all three modes.
- Full gside GTDB zip (3.7 GB) on *A. baumannii*: legacy and gseek JSON are
  byte-identical (`species=Acinetobacter_baumannii`, high,
  `f_unique_weighted=0.5978601252609603`, `ani=98.3`, same two-member
  partition). gseek explicitly passes `threshold_bp=50000` to match the
  legacy sourmash CLI default; the initial 50 bp default admitted extra
  low-weight partitions.
- The benchmark's **rocksdb directory** (`gtdb-k31.rocksdb`) is accepted by
  gseek 0.10.0 as an external reference and hashed as a deterministic tree.
  legacy, gseek, and shadow all returned the same full partition
  (`Acinetobacter_baumannii`, high, `0.5978601252609603`, ANI 98.3); shadow
  produced no difference warning.


## Full-dataset benchmark (`~/temps/test`)

The actual dataset has **42 labels and 4,178 genomes** (the directory includes
the `Shigella` genus label). Every genome was run through legacy and gseek for
all four gseek-backed modes, using 8 worker processes:

| Mode | Runs | Exact same | Errors | Legacy median | gseek median | Speed ratio (median) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| marker | 4,178 | 4,178 | 0 | 0.4381s | 0.4480s | 0.978x |
| panel | 4,178 | 4,178 | 0 | 0.0660s | 0.0748s | 0.882x |
| mash_refseq | 4,178 | 4,178 | 0 | 0.5183s | 0.7013s | 0.739x |
| sourmash | 4,178 | 4,178 | 0 | 1.6415s | 1.6452s | 0.998x |

Total measured time: legacy 11,309s, gseek 12,183s (gseek ≈ 1.08x slower
under the parallel benchmark). The 0.6.0 run reuses ToolResolver probes,
caches reference digests per process, skips manifest/result hashing when not
needed, and uses mash's direct FASTA query path. The benchmark is wall-clock
under 8-way load, not an isolated single-core comparison.

Raw rows and the generated report:
`~/temps/test/results/gseek-full-benchmark-v4/`.
