"""
Downloader for raw genomic FASTQ data from the ENA (European Nucleotide Archive) API.

Adapted for the Escherichia coli ciprofloxacin AMR cohort (Arcadia Science / BV-BRC):
- Supports pilot cohort (100 isolates: 50 R / 50 S, 101 SRA runs)
- Supports full candidate cohort (5,948 isolates: 1,291 R / 4,657 S, 5,964 SRA runs)
- Properly expands isolates with multiple comma-separated SRA runs (e.g. SRR3112315,SRR3878674)
- Treats Genome.ID strictly as a string (prevents float conversion/precision loss)
- Atomic and resumable downloads via temporary files
- Backward-compatible with legacy TSV metadata
"""

import argparse
import concurrent.futures
import csv
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import requests


@dataclass
class SampleRecord:
    """Represents a single bacterial isolate with its genomic and AMR metadata."""
    genome_name: str
    sra_accessions: List[str]
    ciprofloxacin: str
    label: str
    genome_id: str  # Always stored as string to avoid precision loss or ambiguous joins
    biosample_accession: str = ""
    assembly_accession: str = ""
    genome_quality: str = ""
    contigs: str = ""
    contig_l50: str = ""
    contig_n50: str = ""
    checkm_completeness: str = ""
    checkm_contamination: str = ""
    mlst: str = ""
    host_name: str = ""
    isolation_source: str = ""
    geographic_location: str = ""
    extra_fields: Dict[str, str] = field(default_factory=dict)


@dataclass
class SraMetadata:
    """Container for SRA metadata, isolate records, and run mappings."""
    ids: List[str] = field(default_factory=list)  # List of all individual SRA run IDs to download
    samples: List[SampleRecord] = field(default_factory=list)  # List of isolate sample records
    run_to_sample: Dict[str, SampleRecord] = field(default_factory=dict)  # Maps SRA run ID -> SampleRecord

    # Backwards-compatibility fields with legacy pipeline
    mems: List[str] = field(default_factory=list)
    real_mems: List[str] = field(default_factory=list)

    @classmethod
    def from_csv(cls, csv_file: str | Path, delimiter: str = ";") -> "SraMetadata":
        """
        Parses an E. coli metadata CSV file (such as the pilot or candidates CSV).
        Handles UTF-8 with BOM (utf-8-sig), treats Genome.ID as string,
        and expands multiple comma-separated SRA runs per isolate.
        """
        metadata = cls()
        csv_path = Path(csv_file)
        if not csv_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {csv_path}")

        with open(csv_path, mode="r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                clean_row = {
                    k.strip(): (v.strip() if v is not None else "")
                    for k, v in row.items()
                    if k is not None
                }

                genome_name = clean_row.get("Genome.Name", "")
                sra_raw = clean_row.get("SRA.Accession", "")
                cipro = clean_row.get("ciprofloxacin", "")
                label = clean_row.get("label", "")
                # Critical: Genome ID must always remain a string (do not cast to float/int)
                genome_id = str(clean_row.get("Genome.ID", "")).strip()

                # Split multiple SRA runs if present (e.g. "SRR3112315,SRR3878674")
                runs = [
                    run.strip()
                    for run in sra_raw.replace(";", ",").split(",")
                    if run.strip()
                ]

                sample = SampleRecord(
                    genome_name=genome_name,
                    sra_accessions=runs,
                    ciprofloxacin=cipro,
                    label=label,
                    genome_id=genome_id,
                    biosample_accession=clean_row.get("BioSample.Accession", ""),
                    assembly_accession=clean_row.get("Assembly.Accession", ""),
                    genome_quality=clean_row.get("Genome.Quality", ""),
                    contigs=clean_row.get("Contigs", ""),
                    contig_l50=clean_row.get("Contig.L50", ""),
                    contig_n50=clean_row.get("Contig.N50", ""),
                    checkm_completeness=clean_row.get("CheckM.Completeness", ""),
                    checkm_contamination=clean_row.get("CheckM.Contamination", ""),
                    mlst=clean_row.get("MLST", ""),
                    host_name=clean_row.get("Host.Name", ""),
                    isolation_source=clean_row.get("Isolation.Source", ""),
                    geographic_location=clean_row.get("Geographic.Location", ""),
                    extra_fields=clean_row,
                )
                metadata.samples.append(sample)

                for run_id in runs:
                    metadata.ids.append(run_id)
                    metadata.run_to_sample[run_id] = sample
                    metadata.mems.append(label)
                    metadata.real_mems.append(label)

        return metadata

    @classmethod
    def from_tsv(cls, tsv_file: str | Path) -> "SraMetadata":
        """
        Parses a legacy TSV file and initializes the dataclass.
        Skips the header and strips the '.tsv' extension from the ID column.
        """
        metadata = cls()
        with open(tsv_file, mode="r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)

            for row in reader:
                if len(row) >= 3:
                    srr_id = row[0].replace(".tsv", "").strip()
                    metadata.ids.append(srr_id)
                    metadata.mems.append(row[1].strip())
                    metadata.real_mems.append(row[2].strip())
        return metadata

    @classmethod
    def load(cls, file_path: str | Path) -> "SraMetadata":
        """
        Auto-detects file type (.csv or .tsv) and loads metadata.
        """
        path = Path(file_path)
        if path.suffix.lower() == ".tsv":
            return cls.from_tsv(path)
        return cls.from_csv(path)

    def summary(self) -> Dict[str, int]:
        """Returns summary statistics for the loaded metadata."""
        resistant = sum(1 for s in self.samples if s.label.upper() == "R")
        susceptible = sum(1 for s in self.samples if s.label.upper() == "S")
        multi_run_samples = sum(1 for s in self.samples if len(s.sra_accessions) > 1)
        return {
            "total_isolates": len(self.samples),
            "total_runs": len(self.ids),
            "resistant": resistant,
            "susceptible": susceptible,
            "multi_run_isolates": multi_run_samples,
        }

    def export_run_manifest_tsv(self, output_tsv: str | Path) -> None:
        """
        Exports a TSV manifest with one row per SRA run, preserving cross-references
        to Genome Name, Genome ID, BioSample, and resistance phenotype.
        """
        output_path = Path(output_tsv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, mode="w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "ID",
                "Genome_ID",
                "Genome_Name",
                "Ciprofloxacin",
                "Label",
                "BioSample_Accession",
                "Assembly_Accession",
            ])
            for run_id in self.ids:
                sample = self.run_to_sample.get(run_id)
                if sample:
                    writer.writerow([
                        run_id,
                        sample.genome_id,
                        sample.genome_name,
                        sample.ciprofloxacin,
                        sample.label,
                        sample.biosample_accession,
                        sample.assembly_accession,
                    ])
                else:
                    writer.writerow([run_id, "", "", "", "", "", ""])


@dataclass
class FastqDownloader:
    """Handles downloading FASTQ files from the ENA API based on SRA metadata."""

    metadata: SraMetadata
    output_dir: Path = field(default_factory=lambda: Path("database/raw_data"))
    ena_api_url: str = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    max_workers: int = 4
    max_retries: int = 5
    connect_timeout: int = 20
    read_timeout: int = 120
    chunk_size: int = 262144  # 256 KB

    def __post_init__(self):
        self.output_dir = Path(self.output_dir) if self.output_dir else Path("database/raw_data")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def fetch_data(
        self,
        max_bytes: Optional[int] = None,
        limit: Optional[int] = None,
    ) -> None:
        """
        Iterates through the IDs and downloads their corresponding FASTQ files.

        Args:
            max_bytes: If set, only downloads the first N bytes for testing purposes.
            limit: If set, limits download to the first N runs.
        """
        run_ids = self.metadata.ids[:limit] if limit is not None else self.metadata.ids

        if not run_ids:
            logging.info("No SRA runs to download.")
            return

        logging.info(
            f"Starting download for {len(run_ids)} SRA runs ({len(self.metadata.samples)} isolates) "
            f"with {self.max_workers} workers into '{self.output_dir}'..."
        )

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self._process_single_id, run_id, max_bytes): run_id
                for run_id in run_ids
            }

            for future in concurrent.futures.as_completed(futures):
                run_id = futures[future]
                try:
                    future.result()
                except Exception as e:
                    logging.error(f"[{run_id}] A critical error occurred: {e}")

    def _process_single_id(self, run_id: str, max_bytes: Optional[int] = None) -> None:
        sample = self.metadata.run_to_sample.get(run_id)
        sample_info = f" ({sample.genome_name} | {sample.label})" if sample else ""
        logging.debug(f"[{run_id}]{sample_info} Fetching download links...")

        urls = self._get_ena_fastq_urls(run_id)
        if not urls:
            logging.warning(f"[{run_id}] No FASTQ URLs found via ENA API.")
            return

        all_success = True
        for url in urls:
            success = self._download_file(url, max_bytes=max_bytes)
            if not success:
                all_success = False

        if not all_success:
            logging.error(f"[{run_id}] One or more FASTQ files failed to download after all retries.")

    def _get_ena_fastq_urls(self, run_id: str) -> List[str]:
        """
        Queries the ENA API to get direct FTP/HTTP/HTTPS download links for FASTQ files.
        """
        params = {
            "accession": run_id,
            "result": "read_run",
            "fields": "fastq_ftp",
            "format": "json",
        }

        try:
            response = requests.get(self.ena_api_url, params=params, timeout=self.connect_timeout)
            response.raise_for_status()
            data = response.json()

            if not data:
                return []

            ftp_links_str = data[0].get("fastq_ftp", "")
            if not ftp_links_str:
                return []

            urls = []
            for link in ftp_links_str.split(";"):
                link = link.strip()
                if not link:
                    continue
                if link.startswith("ftp://"):
                    link = link[len("ftp://") :]
                if not (link.startswith("http://") or link.startswith("https://")):
                    link = f"https://{link}"
                urls.append(link)
            return urls

        except requests.RequestException as e:
            logging.error(f"Error fetching metadata for {run_id}: {e}")
            return []

    def _download_file(self, url: str, max_bytes: Optional[int] = None) -> bool:
        filename = url.split("/")[-1]

        if max_bytes:
            filename = f"partial_{filename}"

        output_path = self.output_dir / filename

        # If already completely downloaded, skip
        if output_path.exists() and not max_bytes:
            if output_path.stat().st_size > 0:
                logging.debug(f"File {filename} already exists ({output_path.stat().st_size} bytes). Skipping download.")
                return True

        tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")

        # Candidate endpoints: alternate between HTTPS and HTTP if one experiences issues
        urls_to_try = [url]
        if url.startswith("https://ftp.sra.ebi.ac.uk/"):
            urls_to_try.append(url.replace("https://", "http://", 1))
        elif url.startswith("http://ftp.sra.ebi.ac.uk/"):
            urls_to_try.append(url.replace("http://", "https://", 1))

        current_url_idx = 0

        for attempt in range(1, self.max_retries + 1):
            target_url = urls_to_try[current_url_idx % len(urls_to_try)]
            headers = {}
            existing_bytes = 0
            file_mode = "wb"

            # Check if previous partial download exists and resume using Range header
            if tmp_path.exists() and not max_bytes:
                existing_bytes = tmp_path.stat().st_size
                if existing_bytes > 0:
                    headers["Range"] = f"bytes={existing_bytes}-"
                    file_mode = "ab"

            try:
                if existing_bytes > 0:
                    logging.info(
                        f"Downloading {filename} (Attempt {attempt}/{self.max_retries}, resuming from {existing_bytes / (1024*1024):.1f} MB) from {target_url}..."
                    )
                else:
                    logging.info(
                        f"Downloading {filename} (Attempt {attempt}/{self.max_retries}) from {target_url}..."
                    )

                with requests.get(
                    target_url,
                    headers=headers,
                    stream=True,
                    timeout=(self.connect_timeout, self.read_timeout),
                ) as response:
                    # 416: Range Not Satisfiable (e.g. file was already complete or invalid range)
                    if response.status_code == 416:
                        logging.warning(f"[{filename}] Range 416 received. Restarting download from byte 0.")
                        tmp_path.unlink(missing_ok=True)
                        existing_bytes = 0
                        file_mode = "wb"
                        continue

                    # If server ignores Range header (returns 200 instead of 206), overwrite from 0
                    if response.status_code == 200 and existing_bytes > 0:
                        file_mode = "wb"
                        existing_bytes = 0

                    response.raise_for_status()

                    bytes_downloaded = existing_bytes
                    with open(tmp_path, file_mode) as f:
                        for chunk in response.iter_content(chunk_size=self.chunk_size):
                            if chunk:
                                f.write(chunk)
                                bytes_downloaded += len(chunk)

                                if max_bytes and bytes_downloaded >= max_bytes:
                                    logging.debug(f"Reached test limit of {max_bytes} bytes. Stopping download.")
                                    break

                # Atomic finalize
                tmp_path.replace(output_path)
                size_mb = output_path.stat().st_size / (1024 * 1024)
                logging.info(f"Successfully downloaded {filename} ({size_mb:.2f} MB)")
                return True

            except (requests.exceptions.ReadTimeout, requests.exceptions.ConnectionError, requests.exceptions.RequestException) as e:
                logging.warning(
                    f"[{filename}] Connection issue on attempt {attempt}/{self.max_retries}: {e}"
                )
                if len(urls_to_try) > 1:
                    current_url_idx += 1

                if attempt < self.max_retries:
                    wait_time = min(2 * attempt, 30)
                    partial_mb = (tmp_path.stat().st_size / (1024 * 1024)) if tmp_path.exists() else 0
                    logging.info(f"[{filename}] Retrying in {wait_time}s (preserving {partial_mb:.1f} MB already downloaded)...")
                    time.sleep(wait_time)
                else:
                    logging.error(f"[{filename}] Failed after {self.max_retries} attempts.")
                    return False
            except Exception as e:
                logging.error(f"[{filename}] Unexpected error on attempt {attempt}/{self.max_retries}: {e}")
                if attempt < self.max_retries:
                    time.sleep(2)
                else:
                    return False

        return False


def main():
    parser = argparse.ArgumentParser(
        description="Download raw FASTQ data from ENA for E. coli ciprofloxacin AMR cohort."
    )
    parser.add_argument(
        "--file",
        "-f",
        type=str,
        default=None,
        help="Path to specific metadata CSV/TSV file.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download all candidate isolates (5,948 isolates / 5,964 runs) from database/ecoli_ciprofloxacin_RS_candidates_with_metadata.csv.",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="Download the pilot cohort (100 isolates / 101 runs: 50 R, 50 S) from database/ecoli_ciprofloxacin_RS_pilot_50R_50S.csv (DEFAULT).",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=str,
        default="database/raw_data",
        help="Output directory to save downloaded FASTQs (default: database/raw_data).",
    )
    parser.add_argument(
        "--workers",
        "-w",
        type=int,
        default=4,
        help="Number of concurrent download worker threads (default: 4).",
    )
    parser.add_argument(
        "--timeout",
        "-t",
        type=int,
        default=120,
        help="Socket read timeout in seconds for downloading chunks (default: 120s).",
    )
    parser.add_argument(
        "--connect-timeout",
        type=int,
        default=20,
        help="Socket connection timeout in seconds (default: 20s).",
    )
    parser.add_argument(
        "--retries",
        "-r",
        type=int,
        default=5,
        help="Maximum retry attempts on timeout or connection error (default: 5).",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run in test mode (downloads only the first 1MB per file).",
    )
    parser.add_argument(
        "--max-bytes",
        type=int,
        default=None,
        help="Maximum bytes per file to download for testing.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit download to the first N runs.",
    )
    parser.add_argument(
        "--export-manifest",
        type=str,
        default="database/runs_metadata.tsv",
        help="Path to export the expanded run-level TSV manifest (default: database/runs_metadata.tsv).",
    )

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    database_dir = project_root / "database"

    # Select metadata file
    if args.file:
        metadata_path = Path(args.file)
    elif args.all:
        metadata_path = database_dir / "ecoli_ciprofloxacin_RS_candidates_with_metadata.csv"
    else:
        # Default to pilot as recommended in LEIA_PRIMEIRO (Step 3: Baixar os FASTQs do piloto)
        metadata_path = database_dir / "ecoli_ciprofloxacin_RS_pilot_50R_50S.csv"

    if not metadata_path.exists():
        logging.error(f"Metadata file not found: {metadata_path}")
        return

    logging.info(f"Loading metadata from: {metadata_path}")
    metadata = SraMetadata.load(metadata_path)
    summary = metadata.summary()
    logging.info(
        f"Metadata summary: {summary['total_isolates']} isolates "
        f"({summary['resistant']} Resistant, {summary['susceptible']} Susceptible, "
        f"{summary['multi_run_isolates']} with multiple SRA runs) -> "
        f"{summary['total_runs']} total SRA runs to download."
    )

    # Export run-level manifest if requested
    if args.export_manifest:
        manifest_path = (
            project_root / args.export_manifest
            if not Path(args.export_manifest).is_absolute()
            else Path(args.export_manifest)
        )
        metadata.export_run_manifest_tsv(manifest_path)
        logging.info(f"Exported run manifest to {manifest_path}")

    # Set byte limit if in test mode
    max_bytes = args.max_bytes
    if args.test and max_bytes is None:
        max_bytes = 1048576  # 1MB
        logging.info("Test mode enabled: 1MB limit per FASTQ file.")

    out_dir = (
        project_root / args.output_dir
        if not Path(args.output_dir).is_absolute()
        else Path(args.output_dir)
    )

    downloader = FastqDownloader(
        metadata=metadata,
        output_dir=out_dir,
        max_workers=args.workers,
        max_retries=args.retries,
        connect_timeout=args.connect_timeout,
        read_timeout=args.timeout,
    )

    downloader.fetch_data(max_bytes=max_bytes, limit=args.limit)


if __name__ == "__main__":
    main()
