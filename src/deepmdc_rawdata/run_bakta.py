"""
Bakta Annotation Runner & Verifier for DeepMDC Genomic Cohort.
"""

import argparse
import concurrent.futures
import logging
import lzma
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

_pkg_root = Path(__file__).resolve().parent.parent
if str(_pkg_root) not in sys.path:
    sys.path.insert(0, str(_pkg_root))

from deepmdc_rawdata.io_utils import print_summary_banner, read_json_any, resolve_path, write_tsv
from deepmdc_rawdata.metadata import SraMetadata, resolve_metadata_path


@dataclass
class BaktaAnnotationResult:
    """Stores the verification and annotation metrics of a single Bakta run."""
    sample_id: str
    json_path: Path
    exists: bool = False
    is_valid: bool = False
    cds_count: int = 0
    srna_count: int = 0
    trna_count: int = 0
    rrna_count: int = 0
    total_features: int = 0
    error_message: str = ""


class BaktaVerifier:
    """Verifies that Bakta annotation JSON files are complete, valid, and non-empty."""

    @classmethod
    def verify_json(cls, sample_id: str, json_path: Path) -> BaktaAnnotationResult:
        res = BaktaAnnotationResult(sample_id=sample_id, json_path=json_path)
        if not json_path.exists():
            res.error_message = "Annotation JSON file not found."
            return res

        res.exists = True
        if json_path.stat().st_size == 0:
            res.error_message = "Annotation JSON file is empty (0 bytes)."
            return res

        try:
            data = read_json_any(json_path)
            features = data.get("features", [])
            res.total_features = len(features)
            if not features:
                res.error_message = "JSON contains no annotated features."
                return res

            for feat in features:
                ftype = str(feat.get("type", "")).lower()
                if ftype == "cds":
                    res.cds_count += 1
                elif ftype in ("srna", "ncrna"):
                    res.srna_count += 1
                elif ftype == "trna":
                    res.trna_count += 1
                elif ftype == "rrna":
                    res.rrna_count += 1

            if res.cds_count == 0:
                res.error_message = "No CDS features found in annotation."
            else:
                res.is_valid = True
            return res
        except Exception as e:
            res.error_message = f"Error reading/parsing JSON: {e}"
            return res

    @classmethod
    def verify_directory(
        cls,
        sample_ids: List[str],
        jsons_dir: Path,
        report_output: Optional[Path] = None,
    ) -> Tuple[List[BaktaAnnotationResult], List[str]]:
        """Verifies a directory containing Bakta JSONs against a list of sample IDs."""
        results: List[BaktaAnnotationResult] = []
        missing_or_invalid: List[str] = []

        existing_files = {p.name: p for p in jsons_dir.glob("*.*") if p.name.endswith((".json", ".json.xz"))}

        for sid in sample_ids:
            target_path = jsons_dir / f"{sid}.json"
            for cname in (f"{sid}.json.xz", f"{sid}.json"):
                if cname in existing_files:
                    target_path = existing_files[cname]
                    break

            val = cls.verify_json(sid, target_path)
            results.append(val)
            if not val.is_valid:
                missing_or_invalid.append(sid)

        if report_output:
            header = [
                "Sample_ID", "JSON_Path", "Exists", "Is_Valid",
                "CDS_Count", "sRNA_Count", "tRNA_Count", "Total_Features", "Status_Error",
            ]
            rows = [
                [
                    r.sample_id, str(r.json_path), r.exists, r.is_valid,
                    r.cds_count, r.srna_count, r.trna_count, r.total_features, r.error_message or "OK",
                ]
                for r in results
            ]
            write_tsv(report_output, header, rows)
            logging.info(f"Bakta verification report saved to: {report_output}")

        return results, missing_or_invalid


def find_bakta_binary(custom_path: Optional[str] = None) -> Optional[str]:
    """Finds the bakta binary from custom_path, PATH, or common Conda env locations."""
    if custom_path and custom_path != "bakta":
        p = Path(custom_path).expanduser()
        if p.exists() and os.access(p, os.X_OK):
            return str(p.resolve())
        which_custom = shutil.which(custom_path)
        if which_custom:
            return which_custom

    system_bakta = shutil.which("bakta")
    if system_bakta:
        return system_bakta

    home = Path.home()
    conda_prefix = os.environ.get("CONDA_PREFIX")
    candidates = [Path(conda_prefix) / "bin" / "bakta"] if conda_prefix else []
    for folder in ("miniconda3", "anaconda3", ".conda", "miniforge3", "mambaforge"):
        candidates.append(home / folder / "envs" / "bakta" / "bin" / "bakta")

    for cand in candidates:
        if cand.exists() and os.access(cand, os.X_OK):
            return str(cand.resolve())
    return None


@dataclass
class BaktaRunner:
    """Executes Bakta on a set of FASTA genomes with parallelization and auto-compression."""

    fasta_dir: Path
    db_path: Path
    output_dir: Path = field(default_factory=lambda: Path("database/bakta_output"))
    jsons_dir: Path = field(default_factory=lambda: Path("jsons"))
    threads_per_job: int = 8
    parallel_jobs: int = 1
    compress_json: bool = True
    skip_plot: bool = True
    keep_contig_headers: bool = True
    bakta_binary: str = "bakta"

    def __post_init__(self):
        self.fasta_dir = Path(self.fasta_dir)
        self.db_path = Path(self.db_path)
        self.output_dir = Path(self.output_dir)
        self.jsons_dir = Path(self.jsons_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.jsons_dir.mkdir(parents=True, exist_ok=True)

    def is_bakta_available(self) -> bool:
        resolved = find_bakta_binary(self.bakta_binary)
        if resolved:
            self.bakta_binary = resolved
            return True
        return False

    def is_db_valid(self) -> bool:
        return self.db_path.exists() and self.db_path.is_dir()

    def run_single(self, sample_id: str, fasta_path: Path, force: bool = False) -> bool:
        """Runs Bakta on a single genome FASTA file."""
        sample_out_dir = self.output_dir / sample_id
        final_json_path = self.jsons_dir / f"{sample_id}.json"
        final_json_xz_path = self.jsons_dir / f"{sample_id}.json.xz"

        if not force:
            for candidate in (final_json_xz_path, final_json_path):
                if candidate.exists() and BaktaVerifier.verify_json(sample_id, candidate).is_valid:
                    return True

        cmd = [
            self.bakta_binary,
            "--db", str(self.db_path.resolve()),
            "--output", str(sample_out_dir.resolve()),
            "--prefix", sample_id,
            "--threads", str(self.threads_per_job),
            "--force",
        ]
        if self.skip_plot:
            cmd.append("--skip-plot")
        if self.keep_contig_headers:
            cmd.append("--keep-contig-headers")
        cmd.append(str(fasta_path.resolve()))

        logging.info(f"[{sample_id}] Starting Bakta annotation (Threads: {self.threads_per_job})...")
        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if res.returncode != 0:
                logging.error(f"[{sample_id}] Bakta failed (code {res.returncode}):\n{res.stderr[-500:]}")
                return False

            raw_json = sample_out_dir / f"{sample_id}.json"
            if not raw_json.exists():
                logging.error(f"[{sample_id}] Output JSON '{raw_json.name}' was not created.")
                return False

            target = final_json_xz_path if self.compress_json else final_json_path
            if self.compress_json:
                with open(raw_json, "rb") as f_in, lzma.open(target, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            else:
                shutil.copy(raw_json, target)

            val = BaktaVerifier.verify_json(sample_id, target)
            if val.is_valid:
                logging.info(f"[{sample_id}] Annotated: {val.cds_count} CDS, {val.total_features} features.")
                return True
            logging.error(f"[{sample_id}] Verification failed: {val.error_message}")
            return False
        except Exception as e:
            logging.error(f"[{sample_id}] Exception while running Bakta: {e}")
            return False

    def run_cohort(
        self, samples: List[Tuple[str, Path]], force: bool = False
    ) -> Tuple[int, int]:
        """Runs Bakta on a list of (sample_id, fasta_path) tuples."""
        total = len(samples)
        logging.info(
            f"Starting Bakta annotation for {total} genomes ({self.parallel_jobs} job(s), {self.threads_per_job} threads/job)..."
        )
        success_count, fail_count = 0, 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.parallel_jobs) as executor:
            futures = {executor.submit(self.run_single, sid, fpath, force): sid for sid, fpath in samples}
            for idx, future in enumerate(concurrent.futures.as_completed(futures), start=1):
                try:
                    if future.result():
                        success_count += 1
                    else:
                        fail_count += 1
                    logging.info(f"Progress: [{idx}/{total}] (Success: {success_count}, Failed: {fail_count})")
                except Exception as e:
                    logging.error(f"[{futures[future]}] Job crashed: {e}")
                    fail_count += 1
        return success_count, fail_count


def main():
    parser = argparse.ArgumentParser(description="Run and verify Bakta genome annotations for DeepMDC cohort.")
    parser.add_argument("--fasta-dir", "-i", type=str, default="database/fasta_data", help="FASTA directory.")
    parser.add_argument("--db", type=str, default=os.environ.get("BAKTA_DB", ""), help="Bakta database path.")
    parser.add_argument("--bakta-path", "--bakta-bin", type=str, default=os.environ.get("BAKTA_BIN", "bakta"), help="Bakta executable.")
    parser.add_argument("--output-dir", "-o", type=str, default="database/bakta_output", help="Per-genome Bakta output directory.")
    parser.add_argument("--jsons-dir", "-j", type=str, default="jsons", help="Consolidated JSON annotations directory.")
    parser.add_argument("--metadata", "-m", type=str, default=None, help="Path to cohort metadata CSV/TSV.")
    parser.add_argument("--all", action="store_true", help="Use all candidate isolates.")
    parser.add_argument("--pilot", action="store_true", help="Use pilot cohort (100 isolates) (DEFAULT).")
    parser.add_argument("--threads", "-t", type=int, default=8, help="Threads per Bakta execution (default: 8).")
    parser.add_argument("--jobs", type=int, default=1, help="Parallel Bakta jobs (default: 1).")
    parser.add_argument("--no-compress", action="store_true", help="Disable LZMA (.xz) compression.")
    parser.add_argument("--force", action="store_true", help="Re-annotate even if JSON exists.")
    parser.add_argument("--verify", action="store_true", help="Verify existing annotations in --jsons-dir.")
    parser.add_argument("--dry-run", action="store_true", help="List genomes to be annotated without running Bakta.")
    parser.add_argument("--report", type=str, default="database/bakta_verification_report.tsv", help="Verification report TSV.")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    metadata_path = resolve_metadata_path(args.metadata, args.all)
    if not metadata_path.exists():
        logging.error(f"Metadata file not found: {metadata_path}")
        return

    metadata = SraMetadata.load(metadata_path)
    sample_ids = [s.genome_id for s in metadata.samples]
    fasta_dir = resolve_path(args.fasta_dir)
    jsons_dir = resolve_path(args.jsons_dir)
    output_dir = resolve_path(args.output_dir)
    report_path = resolve_path(args.report)

    if args.verify:
        results, missing = BaktaVerifier.verify_directory(sample_ids, jsons_dir, report_output=report_path)
        valid = sum(1 for r in results if r.is_valid)
        total = len(results)
        print_summary_banner(
            "BAKTA ANNOTATION VERIFICATION REPORT",
            [
                ("Total expected isolates", total),
                ("Valid annotations", f"{valid} ({(valid/total*100) if total else 0:.1f}%)"),
                ("Pending or missing", len(missing)),
                ("Detailed report saved to", report_path),
            ],
        )
        return

    if not fasta_dir.exists():
        logging.error(f"FASTA directory '{fasta_dir}' does not exist.")
        return

    existing_fastas = {p.stem: p for p in fasta_dir.glob("*.*") if p.suffix.lower() in (".fna", ".fasta", ".fa")}
    matched_samples: List[Tuple[str, Path]] = []
    for sample in metadata.samples:
        gid = sample.genome_id
        if gid in existing_fastas:
            matched_samples.append((gid, existing_fastas[gid]))
        else:
            for sra in sample.sra_accessions:
                if sra in existing_fastas:
                    matched_samples.append((gid, existing_fastas[sra]))
                    break

    if not matched_samples:
        logging.error(f"No matching FASTA files found in '{fasta_dir}'.")
        return

    if args.dry_run:
        print_summary_banner(
            "DRY-RUN: SAMPLES READY FOR BAKTA ANNOTATION",
            [
                ("Total matched FASTA genomes", len(matched_samples)),
                ("Output directory", output_dir),
                ("Consolidated JSONs directory", jsons_dir),
            ],
        )
        return

    if not args.db:
        logging.error("Bakta database path not provided! Specify --db /path/to/bakta_db or set BAKTA_DB.")
        return

    runner = BaktaRunner(
        fasta_dir=fasta_dir,
        db_path=Path(args.db),
        output_dir=output_dir,
        jsons_dir=jsons_dir,
        threads_per_job=args.threads,
        parallel_jobs=args.jobs,
        compress_json=not args.no_compress,
        bakta_binary=args.bakta_path,
    )

    if not runner.is_bakta_available():
        logging.error("Bakta executable ('bakta') was not found in PATH or Conda environments.")
        return
    if not runner.is_db_valid():
        logging.error(f"Bakta database directory '{args.db}' was not found.")
        return

    runner.run_cohort(matched_samples, force=args.force)
    results, missing = BaktaVerifier.verify_directory(sample_ids, jsons_dir, report_output=report_path)
    valid = sum(1 for r in results if r.is_valid)
    total = len(results)
    print_summary_banner(
        "FINAL BAKTA ANNOTATION SUMMARY",
        [
            ("Total isolates", total),
            ("Completed & valid annotations", f"{valid} ({(valid/total*100) if total else 0:.1f}%)"),
            ("Pending / failed", len(missing)),
            ("Consolidated JSONs", jsons_dir),
            ("Verification report", report_path),
        ],
    )


if __name__ == "__main__":
    main()
