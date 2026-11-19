# deepmdc-ecoli-pipeline

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
conda create -n deepmdc-ecoli-pipeline python=3.11
conda activate deepmdc-ecoli-pipeline
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

Directly retrieves complete assembled bacterial genomes (`.fna`) from BV-BRC and NCBI for the cohort isolates.

### Download Pilot Cohort (100 isolates: 50 Resistant / 50 Susceptible)

```bash
python src/deepmdc_ecoli_pipeline/get_fasta.py
```

### Download All Candidate Isolates (5,948 isolates)

```bash
python src/deepmdc_ecoli_pipeline/get_fasta.py --all
```

### Integrity Verification and Automated Retry

Verify whether all genomes are downloaded and represent valid FASTA files:

```bash
python src/deepmdc_ecoli_pipeline/get_fasta.py --verify
```

Automatically download only missing or corrupted files:

```bash
python src/deepmdc_ecoli_pipeline/get_fasta.py --retry-missing
```

### Additional Options for `get_fasta.py`

- `--output-dir <path>`: Directory to save downloaded FASTA files (default: `database/fasta_data`).
- `--workers <n>`: Number of concurrent download workers (default: `8`).
- `--report <path>`: Path to export detailed validation TSV report (default: `database/fasta_verification_report.tsv`).
- `--manifest <path>`: Path to export isolate-to-FASTA mapping TSV manifest (default: `database/fasta_manifest.tsv`).
- `--file <path>`: Path to a custom metadata CSV or TSV file.
- `--force`: Force re-download of existing files.

---

## 1b. Raw FASTQ Download (`get_fastq.py`)

If you need raw sequencing reads (`.fastq.gz`) from the European Nucleotide Archive (ENA) instead of pre-assembled FASTA contigs, use `src/deepmdc_ecoli_pipeline/get_fastq.py`. It supports resumable downloads via HTTP `Range` headers, automatic HTTPS/HTTP fallback, and GZIP/FASTQ header integrity verification.

### Download Pilot Cohort FASTQs (100 isolates / 101 SRA runs)

```bash
python src/deepmdc_ecoli_pipeline/get_fastq.py
```

### Download All Candidate FASTQs (5,948 isolates / 5,964 SRA runs)

```bash
python src/deepmdc_ecoli_pipeline/get_fastq.py --all
```

### Integrity Verification and Automated Retry

Verify downloaded FASTQ files against the cohort metadata:

```bash
python src/deepmdc_ecoli_pipeline/get_fastq.py --verify
```

Automatically re-download only missing or corrupted SRA runs:

```bash
python src/deepmdc_ecoli_pipeline/get_fastq.py --retry-missing
```

### Additional Options for `get_fastq.py`

- `--output-dir <path>`: Directory to save downloaded FASTQ files (default: `database/raw_data`).
- `--workers <n>`: Number of concurrent download threads (default: `4`).
- `--limit <n>`: Limit download to the first `N` SRA runs.
- `--test`: Test mode — downloads only the first 1 MB of each FASTQ file (`partial_*.fastq.gz`).
- `--report <path>`: Path to export FASTQ verification report (default: `database/fastq_verification_report.tsv`).
- `--manifest <path>`: Path to export FASTQ status manifest (default: `database/fastq_manifest.tsv`).
- `--export-manifest <path>`: Path to export expanded run-level metadata TSV (default: `database/runs_metadata.tsv`).
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
python src/deepmdc_ecoli_pipeline/run_bakta.py \
  --fasta-dir "database/fasta_data" \
  --db "/path/to/bakta_db" \
  --threads 8
```

### Verification and Dry-Run

Verify the status and completeness of generated annotations without re-running Bakta:

```bash
python src/deepmdc_ecoli_pipeline/run_bakta.py --verify --jsons-dir jsons
```

Simulate the annotation execution (Dry-Run):

```bash
python src/deepmdc_ecoli_pipeline/run_bakta.py \
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

Processes Bakta JSON annotation files into the exact dataset structures expected by the original DeepMDC pipeline ([`create_orf_table.py`](https://github.com/CABGenOrg/deepmdc/blob/main/src/utils/create_orf_table.py)). It automatically extracts and separates **ORFs** (protein-coding CDS) and **Non-ORFs** (tRNAs, rRNAs, ncRNAs, sRNAs, CRISPRs, regulatory elements) into dedicated directories:

### Generated Files and Structure

1. **ORFs (Coding Sequences - Amino Acids)**:
   - `database/complete_tables/<ID>.tsv`: Intermediate table with columns `[locus, type, product, sequence, count]`.
   - `database/orfs/<ID>.tsv`: Final DeepRC input table with columns `[orf, templates]`.
   - `database/all_ORFs.json.xz`: Deduplicated global dictionary mapping each amino acid sequence to `{"gene": ..., "product": ..., "ids": ...}` for attention analysis ([`prediction_and_att_weights.py`](https://github.com/CABGenOrg/deepmdc/blob/main/prediction_and_att_weights.py)).

2. **Non-ORFs (Non-Coding Elements - Nucleotides)**:
   - `database/complete_tables_non_orfs/<ID>.tsv`: Intermediate table with columns `[locus, type, product, sequence, count]`.
   - `database/non_orfs/<ID>.tsv`: Final DeepRC input table with columns `[orf, templates]`.
   - `database/all_non_orfs.json.xz`: Deduplicated global dictionary mapping each nucleotide sequence to `{"type": ..., "gene": ..., "product": ..., "ids": ...}`.

3. **Metadata**:
   - `database/metadata.tsv`: DeepMDC-formatted metadata with columns `[ID, CIP, Real_CIP]` (`R`/`S` labels).
   - `database/metadata_non_orfs.tsv`: Dedicated non-ORF metadata matching `database/non_orfs/<ID>.tsv`.

### Run Dataset Generation

```bash
python src/deepmdc_ecoli_pipeline/build_orf_dataset.py
```

### Additional Options for `build_orf_dataset.py`

- `--metadata <path>`: Path to input cohort metadata (default: `database/runs_metadata.tsv`).
- `--jsons-dir <path>`: Directory containing Bakta JSON files (default: `jsons`).
- `--orfs-dir <path>`: Directory to save final DeepRC ORF input tables (default: `database/orfs`).
- `--complete-tables-dir <path>`: Directory to save intermediate complete ORF tables (default: `database/complete_tables`).
- `--all-orfs-output <path>`: Path to save the deduplicated ORF biological dictionary (default: `database/all_ORFs.json.xz`).
- `--non-orfs-dir <path>`: Directory to save final DeepRC non-ORF input tables (default: `database/non_orfs`).
- `--complete-tables-non-orfs-dir <path>`: Directory to save intermediate complete non-ORF tables (default: `database/complete_tables_non_orfs`).
- `--all-non-orfs-output <path>`: Path to save the deduplicated non-ORF dictionary (default: `database/all_non_orfs.json.xz`).
- `--skip-non-orfs`: Skip non-ORF extraction and generate only ORF datasets.
- `--metadata-output <path>`: Path to save the final DeepMDC `metadata.tsv` (default: `database/metadata.tsv`).
- `--metadata-non-orfs-output <path>`: Path to save the non-ORF `metadata_non_orfs.tsv` (default: `database/metadata_non_orfs.tsv`).
- `--antibiotic <name>`: Antibiotic column name (default: `CIP`).
- `--filter-orfs`: Filter out hypothetical proteins and keep only ORF sequences between 40 and 2,000 aa.
- `--use-sra-id`: Use SRA run accessions (e.g. `ERR4035678.tsv`) instead of BV-BRC Genome IDs (`562.100134.tsv`) as filenames.

---

