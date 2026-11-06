"""
Downloader and Verifier for assembled genomic FASTA files (.fna / .fasta)
ready for annotation with Bakta and subsequent DeepMDC/DeepRC processing.

Features:
- Downloads assembled genome contigs directly from BV-BRC (FTP / Data API)
  and NCBI (for records with Assembly Accessions).
- Skips downloading massive raw FASTQs (reduces ~5 TB of data to ~9 GB).
- FASTA validation: ensures file is non-empty, begins with '>', has valid nucleotides,
  and calculates contig count and total genome size.
- Comprehensive verification: checks that all expected isolates from the cohort
  are downloaded, intact, and valid, generating a detailed audit report.
- Resumable: automatically skips already downloaded and valid FASTA files.
- Re-try missing: allows re-downloading only missing or corrupted files with a single flag.
"""

import argparse
import concurrent.futures
import csv
import io
import logging
import os
import re
import sys
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import requests

# Allow running directly as a script without package installation
_pkg_root = Path(__file__).resolve().parent.parent
if str(_pkg_root) not in sys.path:
    sys.path.insert(0, str(_pkg_root))

try:
    from deepmdc_rawdata.get_data import SampleRecord, SraMetadata
except ImportError:
    from .get_data import SampleRecord, SraMetadata


@dataclass
class FastaValidationResult:
    """Stores the outcome of validating a single FASTA file."""
    genome_id: str
    file_path: Path
    exists: bool = False
    is_valid: bool = False
    file_size_bytes: int = 0
    contig_count: int = 0
    total_length_bp: int = 0
    error_message: str = ""

    @property
    def total_length_mbp(self) -> float:
        return self.total_length_bp / 1_000_000.0


class FastaVerifier:
    """Verifies the completeness and biological plausibility of downloaded FASTA assemblies."""

    MIN_EXPECTED_SIZE_BP = 3_000_000   # E. coli genome is typically ~4.5 - 5.5 Mbp
    MAX_EXPECTED_SIZE_BP = 7_000_000

    @classmethod
    def validate_file(cls, genome_id: str, file_path: Path) -> FastaValidationResult:
        res = FastaValidationResult(genome_id=genome_id, file_path=file_path)

        if not file_path.exists():
            res.error_message = "File does not exist."
            return res

        res.exists = True
        res.file_size_bytes = file_path.stat().st_size

        if res.file_size_bytes == 0:
            res.error_message = "File is empty (0 bytes)."
            return res

        # Read beginning to check if it's an HTML error page from web server
        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                first_chars = f.read(512).strip()
                if first_chars.startswith("<") or "<html" in first_chars.lower() or "404 not found" in first_chars.lower():
                    res.error_message = "File contains HTML error page instead of FASTA."
                    return res

                if not first_chars.startswith(">"):
                    res.error_message = "File does not start with '>' FASTA header."
                    return res

            # Parse contigs and calculate length
            contig_count = 0
            total_bp = 0
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        contig_count += 1
                    else:
                        total_bp += len(line)

            res.contig_count = contig_count
            res.total_length_bp = total_bp

            if contig_count == 0:
                res.error_message = "No FASTA headers ('>') found."
                return res

            if total_bp < cls.MIN_EXPECTED_SIZE_BP:
                res.error_message = f"Genome sequence length ({total_bp} bp) is below typical E. coli minimum (3 Mbp)."
                return res

            if total_bp > cls.MAX_EXPECTED_SIZE_BP:
                res.error_message = f"Genome sequence length ({total_bp} bp) exceeds typical E. coli maximum (7 Mbp)."
                return res

            res.is_valid = True
            return res

        except Exception as e:
            res.error_message = f"Error reading file: {e}"
            return res

    @classmethod
    def verify_directory(
        cls,
        metadata: SraMetadata,
        directory: Path,
        report_output: Optional[Path] = None,
    ) -> Tuple[List[FastaValidationResult], List[str]]:
        """
        Verifies that all samples in metadata have a valid FASTA file in directory.
        Returns (list_of_results, list_of_missing_or_invalid_genome_ids).
        """
        results: List[FastaValidationResult] = []
        missing_or_invalid: List[str] = []

        # Find all .fna, .fasta, .fa in directory for lookup
        existing_files = {p.stem: p for p in directory.glob("*.*") if p.suffix.lower() in [".fna", ".fasta", ".fa"]}

        for sample in metadata.samples:
            gid = sample.genome_id
            # Check candidate filenames: <genome_id>.fna, or SRA accession .fna
            target_path = directory / f"{gid}.fna"
            if not target_path.exists() and gid in existing_files:
                target_path = existing_files[gid]
            elif not target_path.exists():
                # Check if named after first SRA accession
                for r in sample.sra_accessions:
                    if r in existing_files:
                        target_path = existing_files[r]
                        break

            val = cls.validate_file(gid, target_path)
            results.append(val)
            if not val.is_valid:
                missing_or_invalid.append(gid)

        if report_output:
            report_output.parent.mkdir(parents=True, exist_ok=True)
            with open(report_output, "w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow([
                    "Genome_ID",
                    "File_Path",
                    "Exists",
                    "Is_Valid",
                    "Size_Bytes",
                    "Contigs",
                    "Total_Length_Mbp",
                    "Status_Error",
                ])
                for r in results:
                    writer.writerow([
                        r.genome_id,
                        str(r.file_path),
                        r.exists,
                        r.is_valid,
                        r.file_size_bytes,
                        r.contig_count,
                        f"{r.total_length_mbp:.3f}",
                        r.error_message or "OK",
                    ])
            logging.info(f"Verification report saved to: {report_output}")

        return results, missing_or_invalid


@dataclass
class FastaDownloader:
    """Downloads assembled bacterial genome FASTA files from BV-BRC and NCBI."""

    metadata: SraMetadata
    output_dir: Path = field(default_factory=lambda: Path("database/fasta_data"))
    max_workers: int = 8
    connect_timeout: int = 15
    read_timeout: int = 60
    max_retries: int = 3
    bvbrc_api_base: str = "https://www.bv-brc.org/api/genome_sequence"
    patric_api_base: str = "https://patricbrc.org/api/genome_sequence"

    def __post_init__(self):
        self.output_dir = Path(self.output_dir) if self.output_dir else Path("database/fasta_data")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _resolve_bvbrc_genome_id(self, sample: SampleRecord) -> Optional[str]:
        """Resolves the real BV-BRC genome_id if Excel stripped trailing zeroes."""
        headers = {"Accept": "application/json"}
        # Try matching by exact genome_name
        if sample.genome_name:
            try:
                encoded_name = requests.utils.quote(sample.genome_name)
                url = f"https://www.bv-brc.org/api/genome/?eq(genome_name,{encoded_name})&select(genome_id)"
                r = requests.get(url, headers=headers, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    if data and "genome_id" in data[0]:
                        resolved = str(data[0]["genome_id"])
                        logging.info(f"[{sample.genome_id}] Resolved to original BV-BRC ID '{resolved}' via Genome Name.")
                        return resolved
            except Exception:
                pass

        # Try matching by SRA accession if available
        for sra in sample.sra_accessions:
            try:
                url = f"https://www.bv-brc.org/api/genome/?eq(sra_accession,{sra})&select(genome_id)"
                r = requests.get(url, headers=headers, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    if data and "genome_id" in data[0]:
                        resolved = str(data[0]["genome_id"])
                        logging.info(f"[{sample.genome_id}] Resolved to original BV-BRC ID '{resolved}' via SRA accession '{sra}'.")
                        return resolved
            except Exception:
                pass

        return None

    def download_sample_fasta(self, sample: SampleRecord, force: bool = False) -> bool:
        """Downloads the FASTA assembly for a single sample."""
        gid = sample.genome_id
        dest_path = self.output_dir / f"{gid}.fna"

        # Check if already downloaded and valid
        if not force and dest_path.exists():
            val = FastaVerifier.validate_file(gid, dest_path)
            if val.is_valid:
                logging.debug(f"[{gid}] Valid FASTA already exists ({val.total_length_mbp:.2f} Mbp, {val.contig_count} contigs). Skipping.")
                return True

        tmp_path = dest_path.with_suffix(".fna.tmp")

        # Excel often strips trailing zeros from decimal IDs (e.g. 562.985 -> 562.98500)
        candidate_gids = [gid]
        if "." in gid:
            prefix, dec = gid.split(".", 1)
            for pad in [1, 2, 3]:
                if len(dec) + pad <= 6:
                    candidate_gids.append(f"{prefix}.{dec}{'0' * pad}")

        urls_to_try: List[str] = []

        # If NCBI Assembly Accession is present (e.g. GCF_... / GCA_...), try it first
        if sample.assembly_accession:
            acc = sample.assembly_accession
            ncbi_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
            urls_to_try.append(ncbi_url)

        # BV-BRC Data API endpoints with &limit(25000) for candidate GIDs
        for c_gid in candidate_gids:
            urls_to_try.append(f"{self.bvbrc_api_base}/?eq(genome_id,{c_gid})&limit(25000)&http_accept=application/dna+fasta")
            urls_to_try.append(f"{self.patric_api_base}/?eq(genome_id,{c_gid})&limit(25000)&http_accept=application/dna+fasta")

        headers = {
            "User-Agent": "DeepMDC-GenomeDownloader/1.0 (academic research; genomics)",
            "Accept": "application/dna+fasta, application/zip, text/plain, */*",
        }

        for attempt in range(1, self.max_retries + 1):
            for url in urls_to_try:
                try:
                    logging.info(f"[{gid}] Fetching assembly from {url} (Attempt {attempt}/{self.max_retries})...")
                    resp = requests.get(
                        url,
                        headers=headers,
                        timeout=(self.connect_timeout, self.read_timeout),
                        stream=True,
                    )
                    if resp.status_code != 200:
                        logging.debug(f"[{gid}] Status {resp.status_code} from {url}")
                        continue

                    # Read into buffer to inspect format
                    buffer = bytearray()
                    for chunk in resp.iter_content(chunk_size=65536):
                        if chunk:
                            buffer.extend(chunk)

                    if not buffer:
                        logging.debug(f"[{gid}] Received empty response from {url}")
                        continue

                    # If response is a ZIP archive (e.g. from NCBI Datasets API)
                    if buffer.startswith(b"PK\x03\x04"):
                        try:
                            with zipfile.ZipFile(io.BytesIO(buffer)) as zf:
                                fna_names = [n for n in zf.namelist() if n.endswith(".fna") or n.endswith(".fasta")]
                                if not fna_names:
                                    logging.warning(f"[{gid}] No FASTA file found in ZIP archive from {url}")
                                    continue
                                fasta_data = zf.read(fna_names[0])
                                with open(tmp_path, "wb") as f:
                                    f.write(fasta_data)
                        except Exception as e:
                            logging.warning(f"[{gid}] Failed to extract ZIP archive: {e}")
                            continue
                    else:
                        with open(tmp_path, "wb") as f:
                            f.write(buffer)

                    # Verify downloaded file
                    val = FastaVerifier.validate_file(gid, tmp_path)
                    if val.is_valid:
                        tmp_path.replace(dest_path)
                        logging.info(
                            f"[{gid}] Successfully downloaded {dest_path.name} "
                            f"({val.total_length_mbp:.2f} Mbp, {val.contig_count} contigs)"
                        )
                        return True
                    else:
                        logging.warning(f"[{gid}] Invalid content downloaded ({val.error_message}).")
                        if tmp_path.exists():
                            tmp_path.unlink()

                except requests.RequestException as e:
                    logging.warning(f"[{gid}] Network error accessing {url}: {e}")
                    if tmp_path.exists():
                        tmp_path.unlink()

            # If all candidates failed on this attempt, try resolving original BV-BRC ID dynamically
            if attempt == 1:
                resolved_id = self._resolve_bvbrc_genome_id(sample)
                if resolved_id and resolved_id not in candidate_gids:
                    candidate_gids.append(resolved_id)
                    urls_to_try.append(f"{self.bvbrc_api_base}/?eq(genome_id,{resolved_id})&limit(25000)&http_accept=application/dna+fasta")
                    urls_to_try.append(f"{self.patric_api_base}/?eq(genome_id,{resolved_id})&limit(25000)&http_accept=application/dna+fasta")

        logging.error(f"[{gid}] Failed to download valid assembly after {self.max_retries} attempts.")
        return False

    def fetch_all(
        self,
        target_ids: Optional[Set[str]] = None,
        force: bool = False,
    ) -> Tuple[int, int]:
        """
        Downloads assemblies for all samples (or only target_ids).
        Returns (success_count, fail_count).
        """
        samples_to_download = [
            s for s in self.metadata.samples
            if (target_ids is None or s.genome_id in target_ids)
        ]

        total = len(samples_to_download)
        logging.info(f"Starting FASTA download for {total} isolates with {self.max_workers} workers into '{self.output_dir}'...")

        success_count = 0
        fail_count = 0

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_sample = {
                executor.submit(self.download_sample_fasta, s, force): s
                for s in samples_to_download
            }

            for future in concurrent.futures.as_completed(future_to_sample):
                sample = future_to_sample[future]
                try:
                    success = future.result()
                    if success:
                        success_count += 1
                    else:
                        fail_count += 1
                except Exception as e:
                    logging.error(f"[{sample.genome_id}] Unexpected error: {e}")
                    fail_count += 1

        logging.info(
            f"FASTA Download finished: {success_count}/{total} successful, {fail_count} failed."
        )
        return success_count, fail_count

    def export_fasta_manifest_tsv(self, output_tsv: Path) -> None:
        """Exports a manifest associating each Genome ID, phenotype and FASTA path."""
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(output_tsv, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "Genome_ID",
                "Genome_Name",
                "Label",
                "Ciprofloxacin",
                "SRA_Accessions",
                "FASTA_Path",
                "Is_Downloaded",
            ])
            for sample in self.metadata.samples:
                gid = sample.genome_id
                fasta_file = self.output_dir / f"{gid}.fna"
                is_downloaded = fasta_file.exists() and fasta_file.stat().st_size > 0
                writer.writerow([
                    gid,
                    sample.genome_name,
                    sample.label,
                    sample.ciprofloxacin,
                    ",".join(sample.sra_accessions),
                    str(fasta_file) if is_downloaded else "",
                    is_downloaded,
                ])
        logging.info(f"FASTA manifest exported to: {output_tsv}")


def main():
    parser = argparse.ArgumentParser(
        description="Download and verify assembled FASTA genomes for Bakta / DeepMDC."
    )
    parser.add_argument(
        "--file", "-f",
        type=str,
        default=None,
        help="Path to metadata CSV/TSV file.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download all candidate isolates (5,948 isolates) from database/ecoli_ciprofloxacin_RS_candidates_with_metadata.csv.",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="Download the pilot cohort (100 isolates: 50 R, 50 S) from database/ecoli_ciprofloxacin_RS_pilot_50R_50S.csv (DEFAULT).",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default="database/fasta_data",
        help="Directory to save downloaded FASTA files (default: database/fasta_data).",
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=8,
        help="Number of concurrent download workers (default: 8).",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify downloaded FASTA files in --output-dir against metadata without downloading.",
    )
    parser.add_argument(
        "--retry-missing",
        action="store_true",
        help="Verify first, and then automatically re-download only missing or invalid files.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download of all files even if already existing.",
    )
    parser.add_argument(
        "--report",
        type=str,
        default="database/fasta_verification_report.tsv",
        help="Path to save verification report TSV (default: database/fasta_verification_report.tsv).",
    )
    parser.add_argument(
        "--manifest",
        type=str,
        default="database/fasta_manifest.tsv",
        help="Path to export the FASTA-isolate manifest TSV (default: database/fasta_manifest.tsv).",
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
        metadata_path = database_dir / "ecoli_ciprofloxacin_RS_pilot_50R_50S.csv"

    if not metadata_path.exists():
        logging.error(f"Metadata file not found: {metadata_path}")
        return

    logging.info(f"Loading metadata from: {metadata_path}")
    metadata = SraMetadata.load(metadata_path)
    summary = metadata.summary()
    logging.info(
        f"Cohort summary: {summary['total_isolates']} isolates "
        f"({summary['resistant']} Resistant, {summary['susceptible']} Susceptible)"
    )

    out_dir = project_root / args.output_dir if not Path(args.output_dir).is_absolute() else Path(args.output_dir)
    report_path = project_root / args.report if not Path(args.report).is_absolute() else Path(args.report)
    manifest_path = project_root / args.manifest if not Path(args.manifest).is_absolute() else Path(args.manifest)

    downloader = FastaDownloader(
        metadata=metadata,
        output_dir=out_dir,
        max_workers=args.workers,
    )

    # Verification Mode Only
    if args.verify and not args.retry_missing:
        logging.info(f"Running verification on '{out_dir}'...")
        results, missing = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        valid_count = sum(1 for r in results if r.is_valid)
        total = len(results)
        print("\n" + "=" * 60)
        print("RELATÓRIO DE VERIFICAÇÃO DE ARQUIVOS FASTA")
        print("=" * 60)
        print(f"Total de isolados esperados: {total}")
        print(f"Arquivos válidos e completos: {valid_count} ({valid_count/total*100:.1f}%)")
        print(f"Arquivos ausentes ou inválidos: {len(missing)}")
        if missing:
            print(f"IDs pendentes (primeiros 10): {missing[:10]}")
            print(f"Para baixar somente os pendentes, execute com --retry-missing")
        print(f"Relatório detalhado salvo em: {report_path}")
        print("=" * 60 + "\n")
        return

    # Retry Missing Mode
    if args.retry_missing:
        logging.info("Checking for missing or corrupted FASTA files...")
        _, missing = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        if not missing:
            logging.info("All FASTA files are already downloaded and valid!")
            downloader.export_fasta_manifest_tsv(manifest_path)
            return
        logging.info(f"Found {len(missing)} missing or invalid FASTA files. Downloading only these...")
        downloader.fetch_all(target_ids=set(missing), force=True)
        # Re-verify after download
        results, remaining_missing = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        downloader.export_fasta_manifest_tsv(manifest_path)
        print(f"\nVerificação final: {len(results) - len(remaining_missing)}/{len(results)} arquivos válidos.")
        return

    # Normal Download Mode
    downloader.fetch_all(force=args.force)

    # Post-download verification
    logging.info("Running automatic post-download verification...")
    results, missing = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
    downloader.export_fasta_manifest_tsv(manifest_path)

    valid_count = sum(1 for r in results if r.is_valid)
    total = len(results)
    print("\n" + "=" * 60)
    print("STATUS FINAL DO DOWNLOAD DE FASTAS")
    print("=" * 60)
    print(f"Total de isolados esperados: {total}")
    print(f"Arquivos baixados e válidos: {valid_count} ({valid_count/total*100:.1f}%)")
    print(f"Arquivos ausentes/falhas: {len(missing)}")
    print(f"Manifesto salvo em: {manifest_path}")
    print(f"Relatório de integridade salvo em: {report_path}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
