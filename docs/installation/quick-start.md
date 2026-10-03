# Quick Start

From zero to species identification in 5 to 15 minutes (including databases).

## Prerequisites

- Linux x86_64
- [pixi](https://pixi.sh) (manages bioinformatics CLI tools)
- Python ≥ 3.12

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

Automatically fetches blast, skani, mash, sourmash, etc.

### 3. Python package

```bash
pip install -e .
```

### 4. Verify

```bash
gside --version
gside db status     # markers should be ✅
```

## Databases

gside uses a tiered database design:

| Tier | Contents | Size | Status |
|------|----------|------|--------|
| mini | Marker rules + sequences | ~1MB | ✅ Bundled with repo |
| panel | Curated reference panel (skani sketch) | ~130MB | Requires `gside db setup` |
| mash | RefSeq MinHash sketch | ~179MB | Requires `gside db setup` |
| all | panel + mash | ~310MB | One-time install |

The `marker` mode works out of the box with the bundled mini tier. The
`panel`, `mash_refseq`, and `all` modes need a one-time database download:

```bash
gside db setup --tier all
```

See [Databases](databases.md) for details.

## Next Steps

- [CLI Reference](../usage/cli.md) — Full command-line options
- [Identification Modes](../usage/modes.md) — Layer principles and use cases
