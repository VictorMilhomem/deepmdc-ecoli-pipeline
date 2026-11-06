"""
Bakta Annotation Runner & Verifier for DeepMDC Genomic Cohort.

Automates running Bakta annotation on assembled bacterial genomes (.fna / .fasta):
- Runs Bakta with optimized parameters (--skip-plot, --keep-contig-headers).
- Manages concurrency: supports multiple parallel jobs and multi-threading per job.
- Resumable: automatically skips already annotated genomes.
- Consolidates JSON outputs into 'jsons/' (or custom directory) ready for
  create_regulatory_elements_table.py.
- Optional compression to .json.xz using LZMA (reduces disk footprint by ~90%,
  fully compatible with handle_regulatory_processing.read_json_any).
- Verifies annotation outputs: counts CDS, sRNAs, and confirms structural integrity.
"""

import argparse
import concurrent.futures
import csv
import json
import logging
import lzma
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Allow running directly as a script
_pkg_root = Path(__file__).resolve().parent.parent
if str(_pkg_root) not in sys.path:
    sys.path.insert(0, str(_pkg_root))

try:
    from deepmdc_rawdata.get_data import SampleRecord, SraMetadata
except ImportError:
    from .get_data import SampleRecord, SraMetadata


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
            if str(json_path).endswith(".xz"):
                with lzma.open(json_path, "rt", encoding="utf-8") as f:
                    data = json.load(f)
            else:
                with open(json_path, "r", encoding="utf-8") as f:
                    data = json.load(f)

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
                return res

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
            # Check sid.json.xz or sid.json
            candidate_names = [f"{sid}.json.xz", f"{sid}.json"]
            target_path = jsons_dir / f"{sid}.json"
            for cname in candidate_names:
                if cname in existing_files:
                    target_path = existing_files[cname]
                    break

            val = cls.verify_json(sid, target_path)
            results.append(val)
            if not val.is_valid:
                missing_or_invalid.append(sid)

        if report_output:
            report_output.parent.mkdir(parents=True, exist_ok=True)
            with open(report_output, "w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow([
                    "Sample_ID",
                    "JSON_Path",
                    "Exists",
                    "Is_Valid",
                    "CDS_Count",
                    "sRNA_Count",
                    "tRNA_Count",
                    "Total_Features",
                    "Status_Error",
                ])
                for r in results:
                    writer.writerow([
                        r.sample_id,
                        str(r.json_path),
                        r.exists,
                        r.is_valid,
                        r.cds_count,
                        r.srna_count,
                        r.trna_count,
                        r.total_features,
                        r.error_message or "OK",
                    ])
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

    # Check system PATH
    system_bakta = shutil.which("bakta")
    if system_bakta:
        return system_bakta

    # Check common Conda environment paths
    home = Path.home()
    conda_prefix = os.environ.get("CONDA_PREFIX")
    conda_candidates = []
    if conda_prefix:
        conda_candidates.append(Path(conda_prefix) / "bin" / "bakta")

    conda_candidates.extend([
        home / "miniconda3" / "envs" / "bakta" / "bin" / "bakta",
        home / "anaconda3" / "envs" / "bakta" / "bin" / "bakta",
        home / ".conda" / "envs" / "bakta" / "bin" / "bakta",
        home / "miniforge3" / "envs" / "bakta" / "bin" / "bakta",
        home / "mambaforge" / "envs" / "bakta" / "bin" / "bakta",
    ])

    for cand in conda_candidates:
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
        """Checks if the bakta binary is available in PATH or Conda."""
        resolved = find_bakta_binary(self.bakta_binary)
        if resolved:
            self.bakta_binary = resolved
            return True
        return False

    def is_db_valid(self) -> bool:
        """Checks if the provided Bakta DB directory is accessible."""
        return self.db_path.exists() and self.db_path.is_dir()

    def run_single(self, sample_id: str, fasta_path: Path, force: bool = False) -> bool:
        """Runs Bakta on a single genome FASTA file."""
        sample_out_dir = self.output_dir / sample_id
        final_json_path = self.jsons_dir / f"{sample_id}.json"
        final_json_xz_path = self.jsons_dir / f"{sample_id}.json.xz"

        # Check if already annotated and valid
        if not force:
            for candidate in [final_json_xz_path, final_json_path]:
                if candidate.exists():
                    val = BaktaVerifier.verify_json(sample_id, candidate)
                    if val.is_valid:
                        logging.debug(f"[{sample_id}] Already annotated ({val.cds_count} CDS, {val.srna_count} sRNA). Skipping.")
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
                logging.error(f"[{sample_id}] Bakta failed with return code {res.returncode}:\n{res.stderr[-500:]}")
                return False

            raw_json = sample_out_dir / f"{sample_id}.json"
            if not raw_json.exists():
                logging.error(f"[{sample_id}] Bakta succeeded but output JSON '{raw_json.name}' was not created.")
                return False

            # Consolidate into jsons_dir
            if self.compress_json:
                logging.info(f"[{sample_id}] Compressing JSON to {final_json_xz_path.name} via LZMA...")
                with open(raw_json, "rb") as f_in, lzma.open(final_json_xz_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                val = BaktaVerifier.verify_json(sample_id, final_json_xz_path)
            else:
                shutil.copy(raw_json, final_json_path)
                val = BaktaVerifier.verify_json(sample_id, final_json_path)

            if val.is_valid:
                logging.info(
                    f"[{sample_id}] Successfully annotated: {val.cds_count} CDS, "
                    f"{val.srna_count} sRNAs, {val.total_features} total features."
                )
                return True
            else:
                logging.error(f"[{sample_id}] Annotation output verification failed: {val.error_message}")
                return False

        except Exception as e:
            logging.error(f"[{sample_id}] Exception while running Bakta: {e}")
            return False

    def run_cohort(
        self,
        samples: List[Tuple[str, Path]],
        force: bool = False,
    ) -> Tuple[int, int]:
        """Runs Bakta on a list of (sample_id, fasta_path) tuples."""
        total = len(samples)
        logging.info(f"Starting Bakta annotation for {total} genomes ({self.parallel_jobs} concurrent job(s), {self.threads_per_job} threads/job)...")

        success_count = 0
        fail_count = 0

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.parallel_jobs) as executor:
            future_to_sample = {
                executor.submit(self.run_single, sid, fpath, force): sid
                for sid, fpath in samples
            }

            for idx, future in enumerate(concurrent.futures.as_completed(future_to_sample), start=1):
                sid = future_to_sample[future]
                try:
                    success = future.result()
                    if success:
                        success_count += 1
                    else:
                        fail_count += 1
                    logging.info(f"Progress: [{idx}/{total}] completed (Success: {success_count}, Failed: {fail_count})")
                except Exception as e:
                    logging.error(f"[{sid}] Job crashed with exception: {e}")
                    fail_count += 1

        return success_count, fail_count


def main():
    parser = argparse.ArgumentParser(
        description="Run and verify Bakta genome annotations for DeepMDC cohort."
    )
    parser.add_argument(
        "--fasta-dir", "-i",
        type=str,
        default="database/fasta_data",
        help="Directory containing assembled FASTA files (.fna / .fasta) (default: database/fasta_data).",
    )
    parser.add_argument(
        "--db",
        type=str,
        default=os.environ.get("BAKTA_DB", ""),
        help="Path to Bakta database directory (or set BAKTA_DB environment variable).",
    )
    parser.add_argument(
        "--bakta-path", "--bakta-bin",
        type=str,
        default=os.environ.get("BAKTA_BIN", "bakta"),
        help="Path to the bakta executable (default: 'bakta' from PATH or auto-detected Conda).",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default="database/bakta_output",
        help="Directory to store detailed per-genome Bakta outputs (default: database/bakta_output).",
    )
    parser.add_argument(
        "--jsons-dir", "-j",
        type=str,
        default="jsons",
        help="Consolidated directory for JSON annotations ready for regulatory table generation (default: jsons).",
    )
    parser.add_argument(
        "--metadata", "-m",
        type=str,
        default="database/ecoli_ciprofloxacin_RS_pilot_50R_50S.csv",
        help="Path to cohort metadata CSV/TSV (default: database/ecoli_ciprofloxacin_RS_pilot_50R_50S.csv).",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Use all candidate isolates from database/ecoli_ciprofloxacin_RS_candidates_with_metadata.csv.",
    )
    parser.add_argument(
        "--pilot",
        action="store_true",
        help="Use pilot cohort (100 isolates) (DEFAULT).",
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=8,
        help="Number of threads per Bakta execution (default: 8).",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of genomes to annotate in parallel (default: 1).",
    )
    parser.add_argument(
        "--no-compress",
        action="store_true",
        help="Do not compress JSON files with LZMA (.xz).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-annotate genomes even if output JSON already exists.",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify existing annotations in --jsons-dir against metadata without running Bakta.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List all genomes and FASTA files that would be annotated without running Bakta.",
    )
    parser.add_argument(
        "--report",
        type=str,
        default="database/bakta_verification_report.tsv",
        help="Path to save verification report TSV (default: database/bakta_verification_report.tsv).",
    )

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent
    database_dir = project_root / "database"

    # Select metadata path
    if args.all:
        metadata_path = database_dir / "ecoli_ciprofloxacin_RS_candidates_with_metadata.csv"
    else:
        metadata_path = Path(args.metadata) if Path(args.metadata).is_absolute() else project_root / args.metadata

    if not metadata_path.exists():
        logging.error(f"Metadata file not found: {metadata_path}")
        return

    metadata = SraMetadata.load(metadata_path)
    sample_ids = [s.genome_id for s in metadata.samples]
    logging.info(f"Loaded {len(sample_ids)} isolates from metadata: {metadata_path}")

    fasta_dir = Path(args.fasta_dir) if Path(args.fasta_dir).is_absolute() else project_root / args.fasta_dir
    jsons_dir = Path(args.jsons_dir) if Path(args.jsons_dir).is_absolute() else project_root / args.jsons_dir
    output_dir = Path(args.output_dir) if Path(args.output_dir).is_absolute() else project_root / args.output_dir
    report_path = Path(args.report) if Path(args.report).is_absolute() else project_root / args.report

    # Verification Mode Only
    if args.verify:
        logging.info(f"Running verification on annotations in '{jsons_dir}'...")
        results, missing = BaktaVerifier.verify_directory(sample_ids, jsons_dir, report_output=report_path)
        valid_count = sum(1 for r in results if r.is_valid)
        total = len(results)
        print("\n" + "=" * 60)
        print("RELATÓRIO DE VERIFICAÇÃO DE ANOTAÇÕES BAKTA")
        print("=" * 60)
        print(f"Total de isolados esperados: {total}")
        print(f"Anotações válidas e completas: {valid_count} ({valid_count/total*100:.1f}%)")
        print(f"Anotações pendentes ou ausentes: {len(missing)}")
        if missing:
            print(f"IDs pendentes (primeiros 10): {missing[:10]}")
        print(f"Relatório detalhado salvo em: {report_path}")
        print("=" * 60 + "\n")
        return

    # Check for FASTA files
    if not fasta_dir.exists():
        logging.error(f"FASTA directory '{fasta_dir}' does not exist. Please download FASTAs first.")
        return

    existing_fastas = {p.stem: p for p in fasta_dir.glob("*.*") if p.suffix.lower() in [".fna", ".fasta", ".fa"]}

    matched_samples: List[Tuple[str, Path]] = []
    missing_fastas: List[str] = []

    for sample in metadata.samples:
        gid = sample.genome_id
        if gid in existing_fastas:
            matched_samples.append((gid, existing_fastas[gid]))
        else:
            # Check if named by SRA accession
            found = False
            for sra in sample.sra_accessions:
                if sra in existing_fastas:
                    matched_samples.append((gid, existing_fastas[sra]))
                    found = True
                    break
            if not found:
                missing_fastas.append(gid)

    if missing_fastas:
        logging.warning(f"{len(missing_fastas)} isolates do not have a FASTA file in '{fasta_dir}' (e.g. {missing_fastas[:5]}).")

    if not matched_samples:
        logging.error(f"No matching FASTA files found in '{fasta_dir}'.")
        return

    logging.info(f"Ready to process {len(matched_samples)} FASTA genomes.")

    if args.dry_run:
        print("\n" + "=" * 60)
        print("DRY-RUN: AMOSTRAS QUE SERÃO ANOTADAS PELO BAKTA")
        print("=" * 60)
        for sid, fpath in matched_samples[:10]:
            print(f"- Sample ID: {sid} | FASTA: {fpath}")
        if len(matched_samples) > 10:
            print(f"... e mais {len(matched_samples) - 10} amostras.")
        print(f"Total: {len(matched_samples)} amostras.")
        print(f"Diretório de saída: {output_dir}")
        print(f"Diretório de JSONs consolidados: {jsons_dir}")
        print("=" * 60 + "\n")
        return

    # Pre-flight check on Bakta binary and database
    if not args.db:
        logging.error(
            "Bakta database path not provided! Please specify --db /path/to/bakta_db "
            "or export BAKTA_DB=/path/to/bakta_db.\n"
            "To download the database, run:\n"
            "  bakta_db download --type light --output <pasta_db>  # (versão leve, ~3.5 GB)\n"
            "  bakta_db download --type full --output <pasta_db>   # (versão completa, ~70 GB)"
        )
        return

    runner = BaktaRunner(
        fasta_dir=fasta_dir,
        db_path=Path(args.db),
        output_dir=output_dir,
        jsons_dir=jsons_dir,
        threads_per_job=args.threads,
        parallel_jobs=args.jobs,
        compress_json=not args.no_compress,
        skip_plot=True,
        keep_contig_headers=True,
        bakta_binary=args.bakta_path,
    )

    if not runner.is_bakta_available():
        logging.error(
            "Bakta executable ('bakta') was not found in PATH!\n"
            "Please install Bakta via Conda/Mamba:\n"
            "  mamba create -n bakta -c conda-forge -c bioconda bakta\n"
            "  conda activate bakta\n"
            "Or run via Docker: oschwengers/bakta"
        )
        return

    if not runner.is_db_valid():
        logging.error(f"Bakta database directory '{args.db}' was not found or is not a directory.")
        return

    success, fail = runner.run_cohort(matched_samples, force=args.force)

    # Post-run verification
    logging.info("Running post-annotation verification...")
    results, missing = BaktaVerifier.verify_directory(sample_ids, jsons_dir, report_output=report_path)
    valid_count = sum(1 for r in results if r.is_valid)
    total = len(results)

    print("\n" + "=" * 60)
    print("STATUS FINAL DA ANOTAÇÃO BAKTA")
    print("=" * 60)
    print(f"Total de isolados: {total}")
    print(f"Anotações concluídas e válidas: {valid_count} ({valid_count/total*100:.1f}%)")
    print(f"Anotações pendentes/falhas: {len(missing)}")
    print(f"JSONs consolidados em: {jsons_dir}")
    print(f"Relatório de anotação salvo em: {report_path}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()

