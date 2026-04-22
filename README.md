# bioinfo-util

A small collection of bioinformatics command-line utilities written in C++.

## Tools

### fastq_stats

Computes per-file statistics from FASTQ or gzipped FASTQ input. Reports yields, read counts, GC content, Q20/Q30 ratios, length distribution (min/avg/median/max), N50, and coverage. Supports optional per-cutoff coverage columns for long-read datasets. Outputs a tab-delimited table to stdout, with timing and memory usage to stderr.

See [`docs/fastq_stats.md`](docs/fastq_stats.md) for full documentation.

## Build and install

### Requirements

- C++17 compiler (GCC 7+ or Clang 5+)
- zlib (`zlib1g-dev` on Debian/Ubuntu, `zlib-devel` on RHEL/Fedora)
- Autotools: `autoconf`, `automake`

### Steps

```bash
autoreconf -i
./configure --prefix=/usr/local
make
make install
```

To install to a custom prefix:

```bash
./configure --prefix=$HOME/.local
make
make install
```

The binary is installed to `$prefix/bin/`. To build without installing:

```bash
autoreconf -i
./configure
make
./src/fastq_stats --help
```
