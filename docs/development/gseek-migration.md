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

## Migrated paths

| Path | Legacy | gseek |
| --- | --- | --- |
| Marker BLAST | `BlastBackend` + `markers_blastdb` | `gseek.align(backend="blast", target=markers.fasta)` |
| FASTQ marker screening | minimap2 PAF + local aggregation | `gseek.align(backend="minimap2", ...)` + same aggregation |

The gseek path deliberately preserves gside semantics:

- coverage is still `aligned_length / length` (not gseek span coverage);
- screening identity is still `nmatch / alen`, via `record.matches`;
- BLAST uses `task="megablast"`, `word_size=11`, and dust on to match the
  historical `blastn` defaults;
- tool discovery uses the same `GSIDE_PIXI_BIN`/pixi/PATH order.

## Not migrated yet

`panel` (skani), `mash_refseq`, and `sourmash` still use their prebuilt database
formats (skani database directories, mash `.msh`, GTDB `.zip`). gseek's index
model currently consumes FASTA and builds its own cache, so switching those
paths needs an explicit external-index API or a gside-side re-sketch step.

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

`shadow` produced no difference warnings. Unit suite: 724 passed; integration:
22 passed, 1 skipped; gseek parity: 2 passed; coverage 97.64%.

