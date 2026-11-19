"""
Downloader and Verifier for assembled genomic FASTA files (.fna / .fasta)
ready for annotation with Bakta and subsequent DeepMDC/DeepRC processing.
"""

import argparse
import concurrent.futures
import io
import logging
import sys
import zipfile
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

        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                first_chars = f.read(512).strip()
                if first_chars.startswith("<") or "<html" in first_chars.lower() or "404 not found" in first_chars.lower():
                    res.error_message = "File contains HTML error page instead of FASTA."
                    return res
                if not first_chars.startswith(">"):
                    res.error_message = "File does not start with '>' FASTA header."
                    return res

            contig_count, total_bp = 0, 0
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
            elif total_bp < cls.MIN_EXPECTED_SIZE_BP:
                res.error_message = f"Genome sequence length ({total_bp} bp) is below typical E. coli minimum (3 Mbp)."
            elif total_bp > cls.MAX_EXPECTED_SIZE_BP:
                res.error_message = f"Genome sequence length ({total_bp} bp) exceeds typical E. coli maximum (7 Mbp)."
            else:
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
        """Verifies that all samples in metadata have a valid FASTA file in directory."""
        results: List[FastaValidationResult] = []
        missing_or_invalid: List[str] = []

        existing_files = {
            p.stem: p for p in directory.glob("*.*") if p.suffix.lower() in (".fna", ".fasta", ".fa")
        }

        for sample in metadata.samples:
            gid = sample.genome_id
            target_path = directory / f"{gid}.fna"
            if not target_path.exists():
                if gid in existing_files:
                    target_path = existing_files[gid]
                else:
                    for r in sample.sra_accessions:
                        if r in existing_files:
                            target_path = existing_files[r]
                            break

            val = cls.validate_file(gid, target_path)
            results.append(val)
            if not val.is_valid:
                missing_or_invalid.append(gid)

        if report_output:
            header = [
                "Genome_ID", "File_Path", "Exists", "Is_Valid",
                "Size_Bytes", "Contigs", "Total_Length_Mbp", "Status_Error",
            ]
            rows = [
                [
                    r.genome_id,
                    str(r.file_path),
                    r.exists,
                    r.is_valid,
                    r.file_size_bytes,
                    r.contig_count,
                    f"{r.total_length_mbp:.3f}",
                    r.error_message or "OK",
                ]
                for r in results
            ]
            write_tsv(report_output, header, rows)
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
        queries = []
        if sample.genome_name:
            queries.append(f"eq(genome_name,{requests.utils.quote(sample.genome_name)})")
        for sra in sample.sra_accessions:
            queries.append(f"eq(sra_accession,{sra})")

        for q in queries:
            try:
                url = f"https://www.bv-brc.org/api/genome/?{q}&select(genome_id)"
                r = requests.get(url, headers=headers, timeout=10)
                if r.status_code == 200:
                    data = r.json()
                    if data and "genome_id" in data[0]:
                        return str(data[0]["genome_id"])
            except Exception:
                pass
        return None

    def download_sample_fasta(self, sample: SampleRecord, force: bool = False) -> bool:
        """Downloads the FASTA assembly for a single sample."""
        gid = sample.genome_id
        dest_path = self.output_dir / f"{gid}.fna"

        if not force and dest_path.exists() and FastaVerifier.validate_file(gid, dest_path).is_valid:
            return True

        tmp_path = dest_path.with_suffix(".fna.tmp")
        candidate_gids = [gid]
        if "." in gid:
            prefix, dec = gid.split(".", 1)
            for pad in (1, 2, 3):
                if len(dec) + pad <= 6:
                    candidate_gids.append(f"{prefix}.{dec}{'0' * pad}")

        urls_to_try: List[str] = []
        if sample.assembly_accession:
            urls_to_try.append(
                f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{sample.assembly_accession}/download?include_annotation_type=GENOME_FASTA"
            )
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
                    resp = requests.get(
                        url, headers=headers, timeout=(self.connect_timeout, self.read_timeout), stream=True
                    )
                    if resp.status_code != 200:
                        continue

                    buffer = bytearray()
                    for chunk in resp.iter_content(chunk_size=65536):
                        if chunk:
                            buffer.extend(chunk)
                    if not buffer:
                        continue

                    if buffer.startswith(b"PK\x03\x04"):
                        with zipfile.ZipFile(io.BytesIO(buffer)) as zf:
                            fna_names = [n for n in zf.namelist() if n.endswith((".fna", ".fasta"))]
                            if not fna_names:
                                continue
                            tmp_path.write_bytes(zf.read(fna_names[0]))
                    else:
                        tmp_path.write_bytes(buffer)

                    val = FastaVerifier.validate_file(gid, tmp_path)
                    if val.is_valid:
                        tmp_path.replace(dest_path)
                        logging.info(f"[{gid}] Downloaded {dest_path.name} ({val.total_length_mbp:.2f} Mbp)")
                        return True
                    tmp_path.unlink(missing_ok=True)
                except Exception as e:
                    logging.warning(f"[{gid}] Error accessing {url}: {e}")
                    tmp_path.unlink(missing_ok=True)

            if attempt == 1:
                resolved_id = self._resolve_bvbrc_genome_id(sample)
                if resolved_id and resolved_id not in candidate_gids:
                    candidate_gids.append(resolved_id)
                    urls_to_try.append(f"{self.bvbrc_api_base}/?eq(genome_id,{resolved_id})&limit(25000)&http_accept=application/dna+fasta")
                    urls_to_try.append(f"{self.patric_api_base}/?eq(genome_id,{resolved_id})&limit(25000)&http_accept=application/dna+fasta")

        logging.error(f"[{gid}] Failed to download valid assembly after {self.max_retries} attempts.")
        return False

    def fetch_all(
        self, target_ids: Optional[Set[str]] = None, force: bool = False
    ) -> Tuple[int, int]:
        """Downloads assemblies for all samples (or only target_ids). Returns (success_count, fail_count)."""
        samples = [s for s in self.metadata.samples if (target_ids is None or s.genome_id in target_ids)]
        logging.info(f"Starting FASTA download for {len(samples)} isolates into '{self.output_dir}'...")

        success_count, fail_count = 0, 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self.download_sample_fasta, s, force): s for s in samples}
            for future in concurrent.futures.as_completed(futures):
                try:
                    if future.result():
                        success_count += 1
                    else:
                        fail_count += 1
                except Exception as e:
                    logging.error(f"[{futures[future].genome_id}] Unexpected error: {e}")
                    fail_count += 1

        return success_count, fail_count

    def export_fasta_manifest_tsv(self, output_tsv: Path) -> None:
        """Exports a manifest associating each Genome ID, phenotype and FASTA path."""
        header = [
            "Genome_ID", "Genome_Name", "Label", "Ciprofloxacin",
            "SRA_Accessions", "FASTA_Path", "Is_Downloaded",
        ]
        rows = []
        for s in self.metadata.samples:
            f = self.output_dir / f"{s.genome_id}.fna"
            ok = f.exists() and f.stat().st_size > 0
            rows.append([
                s.genome_id,
                s.genome_name,
                s.label,
                s.ciprofloxacin,
                ",".join(s.sra_accessions),
                str(f) if ok else "",
                ok,
            ])
        write_tsv(output_tsv, header, rows)


def main():
    parser = argparse.ArgumentParser(description="Download and verify assembled FASTA genomes for Bakta / DeepMDC.")
    parser.add_argument("--file", "-f", type=str, default=None, help="Path to metadata CSV/TSV file.")
    parser.add_argument("--all", action="store_true", help="Download all candidate isolates (5,948 isolates).")
    parser.add_argument("--pilot", action="store_true", help="Download the pilot cohort (100 isolates) (DEFAULT).")
    parser.add_argument("--output-dir", "-o", type=str, default="database/fasta_data", help="Directory for FASTA files.")
    parser.add_argument("--workers", "-w", type=int, default=8, help="Concurrent download workers (default: 8).")
    parser.add_argument("--verify", action="store_true", help="Verify downloaded FASTA files without downloading.")
    parser.add_argument("--retry-missing", action="store_true", help="Verify and re-download only missing/invalid files.")
    parser.add_argument("--force", action="store_true", help="Force re-download of all files.")
    parser.add_argument("--report", type=str, default="database/fasta_verification_report.tsv", help="Verification report TSV.")
    parser.add_argument("--manifest", type=str, default="database/fasta_manifest.tsv", help="FASTA manifest TSV.")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    metadata_path = resolve_metadata_path(args.file, args.all)
    if not metadata_path.exists():
        logging.error(f"Metadata file not found: {metadata_path}")
        return

    metadata = SraMetadata.load(metadata_path)
    out_dir = resolve_path(args.output_dir)
    report_path = resolve_path(args.report)
    manifest_path = resolve_path(args.manifest)

    downloader = FastaDownloader(metadata=metadata, output_dir=out_dir, max_workers=args.workers)

    if args.verify and not args.retry_missing:
        results, missing = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        valid = sum(1 for r in results if r.is_valid)
        total = len(results)
        print_summary_banner(
            "FASTA VERIFICATION REPORT",
            [
                ("Total expected isolates", total),
                ("Valid & complete files", f"{valid} ({(valid/total*100) if total else 0:.1f}%)"),
                ("Missing or invalid files", len(missing)),
                ("Detailed report saved to", report_path),
            ],
        )
        return

    if args.retry_missing:
        _, missing = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        if not missing:
            logging.info("All FASTA files are already downloaded and valid!")
            downloader.export_fasta_manifest_tsv(manifest_path)
            return
        downloader.fetch_all(target_ids=set(missing), force=True)
        results, rem = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
        downloader.export_fasta_manifest_tsv(manifest_path)
        print(f"\nFinal verification: {len(results) - len(rem)}/{len(results)} valid FASTA files.")
        return

    downloader.fetch_all(force=args.force)
    results, missing = FastaVerifier.verify_directory(metadata, out_dir, report_output=report_path)
    downloader.export_fasta_manifest_tsv(manifest_path)
    valid = sum(1 for r in results if r.is_valid)
    total = len(results)
    print_summary_banner(
        "FINAL FASTA DOWNLOAD SUMMARY",
        [
            ("Total expected isolates", total),
            ("Downloaded & valid files", f"{valid} ({(valid/total*100) if total else 0:.1f}%)"),
            ("Missing / failed files", len(missing)),
            ("Manifest saved to", manifest_path),
            ("Integrity report saved to", report_path),
        ],
    )


if __name__ == "__main__":
    main()
