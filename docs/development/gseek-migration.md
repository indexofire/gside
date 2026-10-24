# gseek engine migration

`gside` can route its marker BLAST and FASTQ marker-screening paths through
[`gseek`](https://github.com/indexofire/gseek) behind an explicit switch. The
default remains the existing subprocess implementation; gseek is optional.

## Switch

```bash
GSIDE_ENGINE=legacy   # default: existing backends, no gseek import
GSIDE_ENGINE=gseek    # use gseek's BLAST and minimap2 adapters
GSIDE_ENGINE=shadow   # run both, log differences, return the legacy result
```

`GSIDE_ENGINE` is read per call, so a process can be tested with either mode
without restarting. Invalid values fail before any tool runs.

Install gseek for local testing:

```bash
pip install -e ../gseek
# or, without installing:
PYTHONPATH=../gseek/src pixi run -e dev pytest tests/integration/test_gseek_parity.py -q
```

> The mash/sourmash paths require gseek 0.2.0 (`external_reference`), which is
> published on PyPI: `pip install gseek==0.2.0`. The optional extra is pinned
> accordingly; the sibling checkout remains useful for local development.

## Migrated paths

| Path | Legacy | gseek |
| --- | --- | --- |
| Marker BLAST | `BlastBackend` + `markers_blastdb` | `gseek.align(backend="blast", target=markers.fasta)` |
| FASTQ marker screening | minimap2 PAF + local aggregation | `gseek.align(backend="minimap2", ...)` + same aggregation |
| mash_refseq (assembly/reads) | `MashBackend` against `mash.msh` | `gseek.distance(backend="mash", external_reference=True)` |
| sourmash gather | sketch/gather CSV parsing | `gseek.distance(backend="sourmash", external_reference=True)`, one run per query file |

The gseek path deliberately preserves gside semantics:

- coverage is still `aligned_length / length` (not gseek span coverage);
- screening identity is still `nmatch / alen`, via `record.matches`;
- BLAST uses `task="megablast"`, `word_size=11`, and dust on to match the
  historical `blastn` defaults;
- tool discovery uses the same `GSIDE_PIXI_BIN`/pixi/PATH order.

## Not migrated yet

- `panel` still uses the skani database directory (`panel.sketch`). gseek's
  external-reference contract accepts pre-built files, not a database directory.
- sourmash **rocksdb directories** (`gtdb-k31.rocksdb`, used by the benchmark
  harness for speed) are not accepted by gseek external references. Use the
  GTDB zip/signature file, or keep `GSIDE_ENGINE=legacy`/`shadow` for that
  override. gseek mode now reports an explicit error instead of silently
  returning an empty partition.

## Verification

- `tests/unit/test_gseek_bridge_units.py`: record mapping, aggregation, mode
  dispatch, and fake-gseek parameter wiring.
- `tests/integration/test_gseek_parity.py`: real marker DB/genome comparison
  against the legacy path, plus FASTQ screening parity. Skipped when gseek is
  not importable.

The migration branch must pass the usual gside gate (`ruff`, `pyright`,
`pytest tests/unit --cov ...`). The parity test is an additional integration
check and is not part of the default CI until gseek is an optional dependency
there.


## Local verification (2026-10-24)

`gside scan -m marker -f json` was run against the tracked N16961 genome and a
synthetic FASTQ generated from the `inva` marker:

| Input | legacy | gseek | shadow |
| --- | --- | --- | --- |
| `GCF_000006745.1` assembly | `Vibrio_cholerae` / high / `ompw` 94.767 / 100.1 | identical | identical |
| 21 reads covering `inva` | `Salmonella_enterica` / high / `inva` 100.0 / 98.3 | identical | identical |
| K12 assembly, `-m mash_refseq` | `Escherichia_coli` / high / identity 1.0 | identical | identical |

`shadow` produced no difference warnings. Unit suite: 737 passed; integration:
22 passed, 1 skipped; gseek parity: 5 passed; coverage 97.19%.



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
- A mini sourmash signature DB built from three RefSeq genomes produced
  identical `species`, `confidence`, `f_unique_weighted`, `ani`, and
  `gather_partition` for *A. baumannii* and *S. enterica* in all three modes.
- The benchmark's **rocksdb directory** (`gtdb-k31.rocksdb`) is not a valid
  gseek external reference. legacy still returns the full partition; gseek
  reports an explicit error; shadow returns legacy and logs the error.
