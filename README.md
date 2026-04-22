# fastq_stats

A fast C++ tool for computing per-file FASTQ statistics across both short-read and long-read sequencing data.

It supports plain FASTQ and gzipped FASTQ input, processes multiple files in parallel, reports one line per file, and adds a final `total` row when multiple files are provided.

## Features

- Fast streaming FASTQ parsing in C++
- Supports `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz`
- Per-file TSV output
- Final `total` summary row when multiple input files are given
- Works for both short-read and long-read datasets
- Computes:
  - yields
  - reads
  - GC content
  - Q20 / Q30 base ratios
  - min / average / median / max read length
  - N50
  - coverage
  - optional coverage cutoffs such as `coverage_15k+`, `coverage_100k+`
- File-level parallel processing with `-t / --threads`
- Human-readable size parsing for genome size and coverage cutoffs:
  - `15k`
  - `100k`
  - `1m`
  - `3.1G`

## Output columns

Default output is tab-delimited.

Core columns:

- `file`
- `yields`
- `reads`
- `gc_contents`
- `ratio_Q20`
- `ratio_Q30`
- `min_len`
- `avg_len`
- `med_len`
- `max_len`
- `n50`
- `coverage`

Optional cutoff coverage columns are appended after `coverage`, for example:

- `coverage_15k+`
- `coverage_20k+`
- `coverage_100k+`

Here, `coverage_15k+` means coverage contributed by reads with length **greater than or equal to 15,000 bp**.

## Build

Requires:

- C++17 compiler
- zlib

Example build:

```bash
g++ -O3 -std=c++17 -pthread src/fastq_stats.cpp -lz -o fastq_stats
```

## Usage

```bash
./fastq_stats [options] <fastq1> <fastq2> ...
```

### Options

- `-t, --threads INT`  
  Number of files to process in parallel. Default: `2`

- `--phred INT`  
  Quality score offset. Default: `33`  
  Supported values: `33`, `64`

- `-g, --genome SIZE`  
  Expected genome size for coverage calculation. Default: `3.1G`

- `--cov-cutoff SIZE`  
  Add an extra coverage column for reads with length greater than or equal to the given cutoff.  
  Can be specified multiple times.

- `--no-header`  
  Do not print the header row

- `-h, --help`  
  Show help message

## Examples

Basic usage:

```bash
./fastq_stats sample.fastq.gz
```

Multiple files:

```bash
./fastq_stats sample1.fastq.gz sample2.fastq.gz sample3.fastq.gz
```

Use 4 threads:

```bash
./fastq_stats -t 4 sample1.fastq.gz sample2.fastq.gz
```

Use Phred+64:

```bash
./fastq_stats --phred 64 old_data.fastq.gz
```

Specify genome size:

```bash
./fastq_stats --genome 3.1G sample.fastq.gz
./fastq_stats --genome 3100000000 sample.fastq.gz
```

Add long-read coverage cutoffs:

```bash
./fastq_stats \
  --cov-cutoff 15k \
  --cov-cutoff 20k \
  --cov-cutoff 100k \
  --cov-cutoff 1m \
  sample1.fastq.gz sample2.fastq.gz
```

## Notes

- When multiple files are provided, the program reports one row per file and an additional final row named `total`.
- The `total` row is recalculated from pooled raw statistics, not from averaging per-file summaries.
- Median and N50 are computed using a histogram-based length summary to reduce memory usage while preserving exact results.
- Parallelism is applied across files, not within a single gzipped FASTQ stream.

## Intended use

This tool is designed as a lightweight and fast general-purpose FASTQ statistics utility for:

- Illumina short-read data
- PacBio HiFi data
- ONT data
- ultra-long ONT data

The front half of the output is generally more useful for short-read QC, while the trailing length- and coverage-related fields are especially useful for long-read datasets.

