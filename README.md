# DeepMDC-RawData

This repository provides an automated data retrieval, genome assembly download, Bakta annotation, and feature extraction pipeline for the [DeepMDC](https://github.com/CABGenOrg/deepmdc) project, specifically curated for *Escherichia coli* and Ciprofloxacin antimicrobial resistance (AMR) classification.

---

## Prerequisites

- **Conda** (Miniconda or Anaconda) or **Python >= 3.10**
- **Poetry** (optional, for dependency management)
- **Bakta** (for genome annotation, installed via Conda/Mamba)

---

## Installation

### Using Conda

Create and activate a new Conda environment:

```bash
conda create -n deepmdc-rawdata python=3.11
conda activate deepmdc-rawdata
pip install -r requirements.txt
```

### Using Poetry

Install the project dependencies via Poetry:

```bash
poetry install
```

---

## Pipeline Overview

```
1. Download Assembled Genomes (.fna)       ──►  database/fasta_data/
                   │
2. Genomic Annotation with Bakta           ──►  jsons/*.json.xz
                   │
3. DeepMDC Dataset Generation              ──►  database/orfs/*.tsv
   (ORFs, Complete Tables, all_ORFs, etc.)      database/complete_tables/*.tsv
                                                database/all_ORFs.json.xz
                                                database/metadata.tsv
```

---

## 1. Assembled Genome Download (Recommended)

Directly retrieves complete assembled bacterial genomes (`.fna`) from BV-BRC and NCBI for the cohort isolates. Downloading assembled contigs reduces data transfer from ~5 TB (raw FASTQs) to ~9 GB and eliminates the need for thousands of CPU hours of *de novo* assembly.

### Download Pilot Cohort (100 isolates: 50 Resistant / 50 Susceptible)

```bash
python src/deepmdc_rawdata/get_fasta.py
```

### Download All Candidate Isolates (5,948 isolates)

```bash
python src/deepmdc_rawdata/get_fasta.py --all
```

### Integrity Verification and Automated Retry

Verify whether all genomes are downloaded and represent valid FASTA files:

```bash
python src/deepmdc_rawdata/get_fasta.py --verify
```

Automatically download only missing or corrupted files:

```bash
python src/deepmdc_rawdata/get_fasta.py --retry-missing
```

### Additional Options for `get_fasta.py`

- `--output-dir <path>`: Directory to save downloaded FASTA files (default: `database/fasta_data`).
- `--workers <n>`: Number of concurrent download workers (default: `8`).
- `--report <path>`: Path to export detailed validation TSV report (default: `database/fasta_verification_report.tsv`).
- `--manifest <path>`: Path to export isolate-to-FASTA mapping TSV manifest (default: `database/fasta_manifest.tsv`).
- `--file <path>`: Path to a custom metadata CSV or TSV file.
- `--force`: Force re-download of existing files.

---

## 2. Genomic Annotation with Bakta

Automates batch annotation of assembled FASTA genomes using Bakta, producing consolidated JSON annotations (`.json.xz`) ready for downstream feature extraction.

### Bakta Setup

```bash
# 1. Install Bakta via Conda
conda create -n bakta -c conda-forge -c bioconda bakta
conda activate bakta
pip install requests  # required when running run_bakta.py inside the bakta env

# 2. Download Bakta database (Light database ~3.5 GB recommended for fast high-throughput processing):
bakta_db download --type light --output /path/to/bakta_db
# Or full database (~70 GB):
# bakta_db download --type full --output /path/to/bakta_db
```

> **Tip:** You can set environment variables in your shell so you don't have to specify `--db` every time:
> ```bash
> export BAKTA_DB="/path/to/bakta_db"
> ```
> `run_bakta.py` automatically discovers the `bakta` binary in your PATH or in active/standard Conda environments (`miniconda3/envs/bakta/bin/bakta`).

### Run Cohort Annotation

```bash
python src/deepmdc_rawdata/run_bakta.py \
  --fasta-dir "database/fasta_data" \
  --db "/path/to/bakta_db" \
  --threads 8
```

### Verification and Dry-Run

Verify the status and completeness of generated annotations without re-running Bakta:

```bash
python src/deepmdc_rawdata/run_bakta.py --verify --jsons-dir jsons
```

Simulate the annotation execution (Dry-Run):

```bash
python src/deepmdc_rawdata/run_bakta.py \
  --fasta-dir "database/fasta_data" \
  --dry-run
```

### Additional Options for `run_bakta.py`

- `--threads <n>`: Number of CPU threads per Bakta job (default: `8`).
- `--jobs <n>`: Number of genomes to annotate in parallel (default: `1`).
- `--no-compress`: Disable automatic LZMA compression of output JSONs to `.json.xz`.
- `--output-dir <path>`: Directory to store detailed per-genome Bakta outputs (default: `database/bakta_output`).
- `--jsons-dir <path>`: Consolidated directory for final JSON annotations (default: `jsons`).
- `--force`: Re-annotate genomes even if output JSON already exists.

---

## 3. DeepMDC Dataset Generation (`build_orf_dataset.py`)

Generates the exact dataset structure expected by the original DeepMDC pipeline ([`create_orf_table.py`](https://github.com/CABGenOrg/deepmdc/blob/main/src/utils/create_orf_table.py)):

- `database/complete_tables/<ID>.tsv`: Intermediate table with columns `[locus, type, product, sequence, count]`.
- `database/orfs/<ID>.tsv`: Final DeepRC input table with columns `[orf, templates]`.
- `database/all_ORFs.json.xz`: Deduplicated global dictionary mapping each amino acid sequence to its gene, product, and cross-reference IDs for model attention analysis ([`prediction_and_att_weights.py`](https://github.com/CABGenOrg/deepmdc/blob/main/prediction_and_att_weights.py)).
- `database/metadata.tsv`: DeepMDC-formatted metadata with columns `[ID, CIP, Real_CIP]` (`R`/`S` labels).

### Run Dataset Generation

```bash
python src/deepmdc_rawdata/build_orf_dataset.py
```

### Additional Options for `build_orf_dataset.py`

- `--metadata <path>`: Path to input cohort metadata (default: `database/runs_metadata.tsv`).
- `--jsons-dir <path>`: Directory containing Bakta JSON files (default: `jsons`).
- `--orfs-dir <path>`: Directory to save final DeepRC ORF input tables (default: `database/orfs`).
- `--complete-tables-dir <path>`: Directory to save intermediate complete tables (default: `database/complete_tables`).
- `--all-orfs-output <path>`: Path to save the deduplicated biological dictionary (default: `database/all_ORFs.json.xz`).
- `--metadata-output <path>`: Path to save the final DeepMDC `metadata.tsv` (default: `database/metadata.tsv`).
- `--antibiotic <name>`: Antibiotic column name (default: `CIP`).
- `--filter-orfs`: Filter out hypothetical proteins and keep only sequences between 40 and 2,000 aa.
- `--use-sra-id`: Use SRA run accessions (e.g. `ERR4035678.tsv`) instead of BV-BRC Genome IDs (`562.100134.tsv`) as filenames.

---

## 4. Alternative: Raw FASTQ Download (Unassembled Reads)

If you need the raw unassembled sequencing reads from SRA/ENA:

### Download Pilot Cohort (100 isolates / 101 SRA runs)

```bash
python src/deepmdc_rawdata/get_data.py
```

### Test Mode (1 MB chunk per file)

```bash
python src/deepmdc_rawdata/get_data.py --test
```

### Download All Candidate Isolates (5,948 isolates / 5,964 SRA runs)

```bash
python src/deepmdc_rawdata/get_data.py --all
```

### Additional Options for `get_data.py`

- `--file <path>`: Path to a specific metadata CSV or TSV file.
- `--workers <n>`: Number of concurrent download threads (default: `4`).
- `--output-dir <path>`: Directory to save FASTQ files (default: `database/raw_data`).
- `--timeout <seconds>`: Chunk read timeout in seconds (default: `120s`).
- `--retries <n>`: Maximum retry attempts on timeout or connection failure with HTTP Range resume (default: `5`).
- `--limit <n>`: Limit download to the first N runs.
- `--export-manifest <path>`: Path to export the expanded run-level TSV manifest (default: `database/runs_metadata.tsv`).

---

