# Quick Start

From zero to species identification in 5–15 minutes (including databases).

## Prerequisites

- Linux x86_64
- [pixi](https://pixi.sh) (manages bioinformatics CLI tools)
- Python ≥ 3.11

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/indexofire/gside.git
cd gside
```

### 2. Bioinformatics tools

```bash
pixi install
```

Automatically fetches blast, skani, mash, mmseqs2, etc.

### 3. Python package

```bash
pip install -e .
```

### 4. Verify

```bash
gside --version
gside db status     # markers_v2 should be ✅
```

## Databases

gside uses a tiered database design:

| Tier | Contents | Size | Status |
|------|----------|------|--------|
| mini | Marker rules + sequences | ~1MB | ✅ Bundled with repo |
| panel | Curated reference panel (skani sketch) | ~88MB | Requires `gside db setup` |
| mash | RefSeq MinHash sketch | ~179MB | Requires `gside db setup` |
| all | panel + mash | ~267MB | One-time install |

See [Databases](databases.md) for details.

## Copy from bacmap (local development)

If you already have hermes-bacmap's databases:

```bash
gside db setup --tier all --source /path/to/hermes-bacmap/data/db
```

## Next Steps

- [CLI Reference](../usage/cli.md) — Full command-line options
- [Identification Modes](../usage/modes.md) — Layer principles and use cases
