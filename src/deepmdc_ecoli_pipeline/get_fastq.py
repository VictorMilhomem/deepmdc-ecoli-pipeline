"""
Downloader and Verifier for raw genomic FASTQ files (.fastq.gz) from the
ENA (European Nucleotide Archive) API for Escherichia coli AMR cohorts.
"""

import argparse
import concurrent.futures
import gzip
import logging
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Set, Tuple
import requests

_pkg_root = Path(__file__).resolve().parent.parent
if str(_pkg_root) not in sys.path:
    sys.path.insert(0, str(_pkg_root))

from deepmdc_ecoli_pipeline.io_utils import print_summary_banner, resolve_path, write_tsv
from deepmdc_ecoli_pipeline.metadata import SampleRecord, SraMetadata, resolve_metadata_path


@dataclass
class FastqValidationResult:
    """Stores the validation outcome for a single SRA run's downloaded FASTQ file(s)."""
    run_id: str
    genome_id: str
    file_paths: List[Path] = field(default_factory=list)
    exists: bool = False
    is_valid: bool = False
    file_count: int = 0
    total_size_bytes: int = 0
    error_message: str = ""

    @property
    def total_size_mb(self) -> float:
        return self.total_size_bytes / (1024.0 * 1024.0)


class FastqVerifier:
    """Verifies the presence and GZIP/FASTQ header integrity of downloaded FASTQ files."""

    GZIP_MAGIC = b"\x1f\x8b"

    @classmethod
    def find_run_files(cls, run_id: str, fastq_dir: Path) -> List[Path]:
        """Finds all downloaded FASTQ files matching a given SRA run ID."""
        if not fastq_dir.exists():
            return []
        return sorted(
            p for p in fastq_dir.glob(f"{run_id}*.fastq.gz")
            if p.is_file() and not p.name.endswith(".tmp")
        )

    @classmethod
    def validate_gzip_fastq(cls, file_path: Path) -> Tuple[bool, str]:
        """Validates that a .fastq.gz file is non-empty, has GZIP magic bytes, and starts with '@'."""
        if not file_path.exists():
            return False, "File does not exist."
        if file_path.stat().st_size == 0:
            return False, "File is empty (0 bytes)."

        try:
            with open(file_path, "rb") as f:
                magic = f.read(2)
                if magic != cls.GZIP_MAGIC:
                    return False, f"Invalid GZIP magic header ({magic!r})."

            with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as gz:
                first_line = gz.readline().strip()
                if not first_line:
                    return False, "GZIP stream is empty."
                if not first_line.startswith("@"):
                    return False, "Decompressed content does not start with '@' FASTQ header."
        except (gzip.BadGzipFile, EOFError, OSError) as e:
            return False, f"Corrupted GZIP stream: {e}"

        return True, ""

    @classmethod
    def validate_run(cls, run_id: str, genome_id: str, fastq_dir: Path) -> FastqValidationResult:
        """Validates all FASTQ files belonging to a single SRA run ID."""
        files = cls.find_run_files(run_id, fastq_dir)
        res = FastqValidationResult(
            run_id=run_id,
            genome_id=genome_id,
            file_paths=files,
            file_count=len(files),
        )
        if not files:
            res.error_message = "No FASTQ files found for run."
            return res

        res.exists = True
        res.total_size_bytes = sum(f.stat().st_size for f in files)
        for fpath in files:
            ok, err = cls.validate_gzip_fastq(fpath)
            if not ok:
                res.error_message = f"{fpath.name}: {err}"
                return res

        res.is_valid = True
        return res

    @classmethod
    def verify_directory(
        cls,
        metadata: SraMetadata,
        fastq_dir: str | Path,
        report_output: Optional[str | Path] = None,
    ) -> Tuple[List[FastqValidationResult], List[str]]:
        """Audits all expected SRA runs in metadata against fastq_dir."""
        out_dir = Path(fastq_dir)
        results: List[FastqValidationResult] = []
        missing_or_invalid: List[str] = []

        for run_id in metadata.ids:
            sample = metadata.run_to_sample.get(run_id)
            res = cls.validate_run(run_id, sample.genome_id if sample else "", out_dir)
            results.append(res)
            if not res.is_valid:
                missing_or_invalid.append(run_id)

        if report_output:
            header = [
                "Run_ID", "Genome_ID", "Genome_Name", "Label",
                "Exists", "Is_Valid", "File_Count", "Total_Size_MB", "Files", "Error_Message",
            ]
            rows = []
            for res in results:
                s = metadata.run_to_sample.get(res.run_id)
                rows.append([
                    res.run_id,
                    res.genome_id,
                    s.genome_name if s else "",
                    s.label if s else "",
                    res.exists,
                    res.is_valid,
                    res.file_count,
                    f"{res.total_size_mb:.2f}",
                    ";".join(p.name for p in res.file_paths),
                    res.error_message,
                ])
            write_tsv(report_output, header, rows)

        return results, missing_or_invalid


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
        target_ids: Optional[Set[str]] = None,
        force: bool = False,
    ) -> None:
        """Iterates through SRA run IDs and downloads their corresponding FASTQ files."""
        run_ids = [r for r in self.metadata.ids if (target_ids is None or r in target_ids)]
        if limit is not None:
            run_ids = run_ids[:limit]

        if not run_ids:
            logging.info("No SRA runs to download.")
            return

        logging.info(
            f"Starting download for {len(run_ids)} SRA runs ({len(self.metadata.samples)} isolates) "
            f"with {self.max_workers} workers into '{self.output_dir}'..."
        )

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self._process_single_id, run_id, max_bytes, force): run_id
                for run_id in run_ids
            }
            for future in concurrent.futures.as_completed(futures):
                run_id = futures[future]
                try:
                    future.result()
                except Exception as e:
                    logging.error(f"[{run_id}] A critical error occurred: {e}")

    def _process_single_id(
        self, run_id: str, max_bytes: Optional[int] = None, force: bool = False
    ) -> None:
        urls = self._get_ena_fastq_urls(run_id)
        if not urls:
            logging.warning(f"[{run_id}] No FASTQ URLs found via ENA API.")
            return

        if not all(self._download_file(url, max_bytes=max_bytes, force=force) for url in urls):
            logging.error(f"[{run_id}] One or more FASTQ files failed to download after all retries.")

    def _get_ena_fastq_urls(self, run_id: str) -> List[str]:
        """Queries the ENA API to get direct HTTP/HTTPS download links for FASTQ files."""
        params = {"accession": run_id, "result": "read_run", "fields": "fastq_ftp", "format": "json"}
        try:
            response = requests.get(self.ena_api_url, params=params, timeout=self.connect_timeout)
            response.raise_for_status()
            data = response.json()
            if not data or not data[0].get("fastq_ftp"):
                return []

            urls = []
            for link in data[0]["fastq_ftp"].split(";"):
                link = link.strip()
                if not link:
                    continue
                if link.startswith("ftp://"):
                    link = link[len("ftp://") :]
                if not link.startswith(("http://", "https://")):
                    link = f"https://{link}"
                urls.append(link)
            return urls
        except requests.RequestException as e:
            logging.error(f"Error fetching metadata for {run_id}: {e}")
            return []

    def _download_file(
        self, url: str, max_bytes: Optional[int] = None, force: bool = False
    ) -> bool:
        filename = f"partial_{url.split('/')[-1]}" if max_bytes else url.split("/")[-1]
        output_path = self.output_dir / filename

        if output_path.exists() and not max_bytes and not force and output_path.stat().st_size > 0:
            logging.debug(f"File {filename} already exists. Skipping download.")
            return True

        tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
        if force and tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

        urls_to_try = [url]
        if url.startswith("https://ftp.sra.ebi.ac.uk/"):
            urls_to_try.append(url.replace("https://", "http://", 1))
        elif url.startswith("http://ftp.sra.ebi.ac.uk/"):
            urls_to_try.append(url.replace("http://", "https://", 1))

        current_url_idx = 0
        for attempt in range(1, self.max_retries + 1):
            target_url = urls_to_try[current_url_idx % len(urls_to_try)]
            headers = {}
            existing_bytes = tmp_path.stat().st_size if (tmp_path.exists() and not max_bytes) else 0
            file_mode = "ab" if existing_bytes > 0 else "wb"
            if existing_bytes > 0:
                headers["Range"] = f"bytes={existing_bytes}-"

            try:
                logging.info(
                    f"Downloading {filename} (Attempt {attempt}/{self.max_retries}) from {target_url}..."
                )
                with requests.get(
                    target_url,
                    headers=headers,
                    stream=True,
                    timeout=(self.connect_timeout, self.read_timeout),
                ) as response:
                    if response.status_code == 416:
                        tmp_path.unlink(missing_ok=True)
                        continue
                    if response.status_code == 200 and existing_bytes > 0:
                        file_mode, existing_bytes = "wb", 0

                    response.raise_for_status()
                    bytes_downloaded = existing_bytes
                    with open(tmp_path, file_mode) as f:
                        for chunk in response.iter_content(chunk_size=self.chunk_size):
                            if chunk:
                                f.write(chunk)
                                bytes_downloaded += len(chunk)
                                if max_bytes and bytes_downloaded >= max_bytes:
                                    break

                tmp_path.replace(output_path)
                size_mb = output_path.stat().st_size / (1024 * 1024)
                logging.info(f"Successfully downloaded {filename} ({size_mb:.2f} MB)")
                return True

            except requests.RequestException as e:
                logging.warning(f"[{filename}] Connection issue on attempt {attempt}/{self.max_retries}: {e}")
                current_url_idx += 1
                if attempt < self.max_retries:
                    time.sleep(min(2 * attempt, 30))
                else:
                    return False
            except Exception as e:
                logging.error(f"[{filename}] Unexpected error on attempt {attempt}/{self.max_retries}: {e}")
                if attempt < self.max_retries:
                    time.sleep(2)
                else:
                    return False

        return False

    def export_fastq_manifest_tsv(self, output_tsv: str | Path) -> None:
        """Exports a TSV manifest mapping each SRA run to its isolate metadata and local FASTQ paths."""
        header = [
            "Run_ID", "Genome_ID", "Genome_Name", "Label",
            "Ciprofloxacin", "FASTQ_Files", "Downloaded",
        ]
        rows = []
        for run_id in self.metadata.ids:
            s = self.metadata.run_to_sample.get(run_id)
            files = FastqVerifier.find_run_files(run_id, self.output_dir)
            is_ok = bool(files) and all(FastqVerifier.validate_gzip_fastq(p)[0] for p in files)
            rows.append([
                run_id,
                s.genome_id if s else "",
                s.genome_name if s else "",
                s.label if s else "",
                s.ciprofloxacin if s else "",
                ";".join(str(p) for p in files) if is_ok else "",
                is_ok,
            ])
        write_tsv(output_tsv, header, rows)
        logging.info(f"FASTQ manifest exported to: {output_tsv}")


def main():
    parser = argparse.ArgumentParser(
        description="Download and verify raw E. coli FASTQ files from ENA for DeepMDC."
    )
    parser.add_argument("--file", "-f", type=str, default=None, help="Path to metadata CSV/TSV file.")
    parser.add_argument("--all", action="store_true", help="Download all candidate isolates (5,948 isolates / 5,964 runs).")
    parser.add_argument("--pilot", action="store_true", help="Download the pilot cohort (100 isolates / 101 runs) (DEFAULT).")
    parser.add_argument("--output-dir", "-o", type=str, default="database/raw_data", help="Output directory for FASTQs.")
    parser.add_argument("--workers", "-w", type=int, default=4, help="Concurrent download threads (default: 4).")
    parser.add_argument("--timeout", "-t", type=int, default=120, help="Read timeout in seconds (default: 120s).")
    parser.add_argument("--connect-timeout", type=int, default=20, help="Connect timeout in seconds (default: 20s).")
    parser.add_argument("--retries", "-r", type=int, default=5, help="Maximum retry attempts (default: 5).")
    parser.add_argument("--verify", action="store_true", help="Verify downloaded FASTQ files against metadata.")
    parser.add_argument("--retry-missing", action="store_true", help="Verify and re-download only missing/corrupted FASTQ runs.")
    parser.add_argument("--force", action="store_true", help="Force re-download of all FASTQ files.")
    parser.add_argument("--test", action="store_true", help="Test mode (downloads only first 1MB per file).")
    parser.add_argument("--max-bytes", type=int, default=None, help="Maximum bytes per file to download for testing.")
    parser.add_argument("--limit", type=int, default=None, help="Limit download to the first N runs.")
    parser.add_argument("--export-manifest", type=str, default="database/runs_metadata.tsv", help="Path to export run-level TSV.")
    parser.add_argument("--manifest", type=str, default="database/fastq_manifest.tsv", help="Path to export FASTQ status manifest.")
    parser.add_argument("--report", type=str, default="database/fastq_verification_report.tsv", help="Path to save verification report.")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    metadata_path = resolve_metadata_path(args.file, args.all)
    if not metadata_path.exists():
        logging.error(f"Metadata file not found: {metadata_path}")
        return

    metadata = SraMetadata.load(metadata_path)
    s = metadata.summary()
    logging.info(
        f"Loaded {s['total_isolates']} isolates ({s['resistant']} R, {s['susceptible']} S) -> {s['total_runs']} SRA runs."
    )

    out_dir = resolve_path(args.output_dir)
    report_path = resolve_path(args.report)
    fastq_manifest_path = resolve_path(args.manifest)

    if args.export_manifest:
        run_manifest_path = resolve_path(args.export_manifest)
        metadata.export_run_manifest_tsv(run_manifest_path)

    downloader = FastqDownloader(
        metadata=metadata,
        output_dir=out_dir,
        max_workers=args.workers,
        max_retries=args.retries,
        connect_timeout=args.connect_timeout,
        read_timeout=args.timeout,
    )

    if args.verify and not args.retry_missing:
        results, missing = FastqVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        valid = sum(1 for r in results if r.is_valid)
        total = len(results)
        items = [
            ("Total expected SRA runs", total),
            ("Valid & complete runs", f"{valid} ({(valid/total*100) if total else 0:.1f}%)"),
            ("Missing or invalid runs", len(missing)),
            ("Detailed report saved to", report_path),
        ]
        if missing:
            items.insert(3, ("Pending run IDs (first 10)", missing[:10]))
        print_summary_banner("FASTQ VERIFICATION REPORT", items)
        return

    if args.retry_missing:
        _, missing = FastqVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        if not missing:
            logging.info("All FASTQ runs are already downloaded and valid!")
            downloader.export_fastq_manifest_tsv(fastq_manifest_path)
            return
        downloader.fetch_data(target_ids=set(missing), limit=args.limit, force=True)
        results, rem = FastqVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        downloader.export_fastq_manifest_tsv(fastq_manifest_path)
        print(f"\nFinal verification: {len(results) - len(rem)}/{len(results)} valid FASTQ runs.")
        return

    max_bytes = args.max_bytes or (1048576 if args.test else None)
    downloader.fetch_data(max_bytes=max_bytes, limit=args.limit, force=args.force)

    if not max_bytes:
        results, missing = FastqVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        downloader.export_fastq_manifest_tsv(fastq_manifest_path)
        valid = sum(1 for r in results if r.is_valid)
        total = len(results)
        print_summary_banner(
            "FINAL FASTQ DOWNLOAD SUMMARY",
            [
                ("Total expected SRA runs", total),
                ("Downloaded & valid runs", f"{valid} ({(valid/total*100) if total else 0:.1f}%)"),
                ("Missing / failed runs", len(missing)),
                ("Manifest saved to", fastq_manifest_path),
                ("Integrity report saved to", report_path),
            ],
        )


if __name__ == "__main__":
    main()
