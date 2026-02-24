# methfast (Rust)

`methfast` is a Rust CLI for extracting weighted methylation values from a bed-like methylation file over target BED intervals.

## Features

- Plain text and gzipped methylation input support
- Same core behavior and output format as `methfast` C `v0.3.0`
- Compatible short flags plus modern long flags
- Parallel target processing (`--threads`)

## Build

### macOS

1. Install Rust via [rustup](https://rustup.rs/).
2. Build:

```bash
cargo build --release
```

Binary path:

```bash
./target/release/methfast
```

### Linux

1. Install Rust via [rustup](https://rustup.rs/).
2. Build:

```bash
cargo build --release
```

Binary path:

```bash
./target/release/methfast
```

## Install

Install from the current repo:

```bash
cargo install --path . --locked
```

This installs `methfast` into Cargo's bin directory (usually `$HOME/.cargo/bin`).

## Usage

```bash
methfast <methylation_bed(.gz)> <target_bed> [OPTIONS]
```

### Positional arguments

- `METHYLATION_BED`: bedmethyl-style input (`.bed` or `.bed.gz`)
- `TARGET_BED`: target BED intervals

### Options

- `-f, --fraction-col <INT>`: methylation fraction column (1-based, default `4`)
- `-c, --coverage-col <INT>`: total coverage column (1-based, default `5`)
- `-m, --methylated-col <INT>`: methylated coverage column (1-based)
- `-u, --unmethylated-col <INT>`: unmethylated coverage column (1-based)
- `-o, --output <FILE>`: output file (default: stdout)
- `-t, --threads <INT>`: worker thread count for target processing

## Output format

Tab-separated columns:

1. chrom
2. start
3. end
4. number of overlapping methylation positions
5. summed total coverage over overlaps
6. weighted methylation fraction (4 decimals)

## Development checks

```bash
make ci
```

Equivalent commands:

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
```

## Releases (GitHub Actions)

- CI release workflow: `.github/workflows/release.yml`
- Trigger: push a version tag like `v0.1.0`
- Outputs: GitHub Release assets for:
  - `x86_64-unknown-linux-gnu`
  - `aarch64-apple-darwin`
  - `x86_64-apple-darwin`

Example:

```bash
git tag v0.1.0
git push origin v0.1.0
```

The workflow builds with `--locked`, packages the `methfast` binary plus `README.md` and `LICENSE`, and uploads archives to the GitHub Release page.
