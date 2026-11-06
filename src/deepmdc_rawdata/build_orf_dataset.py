"""
DeepMDC ORF Dataset Generator.

Processes Bakta JSON annotation files (.json or .json.xz) into the exact format
expected by DeepMDC (deepmdc/src/utils):
1. database/complete_tables/<ID>.tsv:
   Intermediate complete table containing [locus, type, product, sequence, count].
2. database/orfs/<ID>.tsv:
   Final DeepRC input table containing [orf, templates].
3. database/all_ORFs.json.xz:
   Deduplicated global biological dictionary mapping every amino acid sequence
   to its gene symbol, product annotation, and database cross-reference IDs.
4. database/metadata.tsv:
   Metadata table linking sample filenames to their antibiotic phenotype label
   (e.g., ID: 562.100134.tsv, CIP: R/S, Real_CIP: R/S).
"""

import argparse
import csv
import json
import logging
import lzma
import os
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable


# Add package root to sys.path
_pkg_root = Path(__file__).resolve().parent.parent
if str(_pkg_root) not in sys.path:
    sys.path.insert(0, str(_pkg_root))

try:
    from deepmdc_rawdata.get_data import SampleRecord, SraMetadata
except ImportError:
    from .get_data import SampleRecord, SraMetadata


@dataclass
class OrfExtractionResult:
    """Stores metrics for a single sample's extracted ORFs."""
    sample_id: str
    total_orfs: int
    unique_orfs: int
    complete_table_path: Path
    orf_table_path: Path


class OrfDatasetBuilder:
    """
    Builds the DeepMDC ORF dataset and the consolidated all_ORFs dictionary.
    """

    def __init__(
        self,
        metadata_file: Union[str, Path],
        jsons_dir: Union[str, Path] = "jsons",
        complete_tables_dir: Union[str, Path] = "database/complete_tables",
        orfs_dir: Union[str, Path] = "database/orfs",
        all_orfs_output: Union[str, Path] = "database/all_ORFs.json.xz",
        metadata_output: Union[str, Path] = "database/metadata.tsv",
        antibiotic_name: str = "CIP",
        filter_orfs: bool = False,
        use_sra_id: bool = False,
        min_aa_len: int = 40,
        max_aa_len: int = 2000,
    ):
        self.metadata_file = Path(metadata_file)
        self.jsons_dir = Path(jsons_dir)
        self.complete_tables_dir = Path(complete_tables_dir)
        self.orfs_dir = Path(orfs_dir)
        self.all_orfs_output = Path(all_orfs_output)
        self.metadata_output = Path(metadata_output)
        self.antibiotic_name = antibiotic_name.upper()
        self.filter_orfs = filter_orfs
        self.use_sra_id = use_sra_id
        self.min_aa_len = min_aa_len
        self.max_aa_len = max_aa_len

        # Ensure output directories exist
        self.complete_tables_dir.mkdir(parents=True, exist_ok=True)
        self.orfs_dir.mkdir(parents=True, exist_ok=True)
        self.all_orfs_output.parent.mkdir(parents=True, exist_ok=True)
        self.metadata_output.parent.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def read_compressed_json(json_path: Union[str, Path]) -> dict:
        """Reads a plain (.json) or LZMA-compressed (.json.xz) JSON file."""
        p = Path(json_path)
        if not p.exists():
            return {}
        if str(p).endswith(".xz"):
            with lzma.open(p, "rt", encoding="utf-8") as inp:
                return json.load(inp)
        else:
            with open(p, "r", encoding="utf-8") as inp:
                return json.load(inp)

    def find_json_file(self, genome_id: str, sra_ids: List[str]) -> Optional[Path]:
        """Locates the Bakta JSON file matching either genome_id or any of its SRA run IDs."""
        candidates = []
        # Prioritize genome_id
        candidates.extend([
            self.jsons_dir / f"{genome_id}.json.xz",
            self.jsons_dir / f"{genome_id}.json",
        ])
        # Check SRA IDs
        for sra in sra_ids:
            candidates.extend([
                self.jsons_dir / f"{sra}.json.xz",
                self.jsons_dir / f"{sra}.json",
            ])

        for c in candidates:
            if c.exists() and c.is_file():
                return c
        return None

    def process_single_json(
        self,
        json_file: Path,
    ) -> Tuple[List[List[Union[str, int]]], Dict[str, dict]]:
        """
        Extracts CDS amino acid sequences from a Bakta JSON file.

        Returns:
            Tuple of:
            - rows: List of [locus, type, product, sequence, count]
            - orfs_meta: Dict mapping aa_sequence -> {"gene": str|None, "product": str, "ids": str}
        """
        data = self.read_compressed_json(json_file)
        features = data.get("features", [])

        sequences: List[str] = []
        raw_rows: List[List[str]] = []
        orfs_meta: Dict[str, dict] = {}

        for feat in features:
            feat_type = str(feat.get("type", "")).lower()
            if feat_type != "cds":
                continue

            sequence = str(feat.get("aa", "")).strip().upper()
            if not sequence:
                continue

            product = str(feat.get("product", "")).strip()
            locus = str(feat.get("locus", feat.get("id", ""))).strip()
            gene = feat.get("gene") or None
            if gene:
                gene = str(gene).strip()

            if self.filter_orfs:
                if product.lower() == "hypothetical protein":
                    continue
                if not (self.min_aa_len <= len(sequence) <= self.max_aa_len):
                    continue

            sequences.append(sequence)
            raw_rows.append([locus, "cds", product, sequence])

            # Deduplicate biological annotation into dictionary
            if sequence not in orfs_meta:
                db_xrefs = feat.get("db_xrefs", [])
                if isinstance(db_xrefs, list):
                    ids_str = ", ".join(str(x) for x in db_xrefs)
                else:
                    ids_str = str(db_xrefs or "")

                orfs_meta[sequence] = {
                    "gene": gene,
                    "product": product,
                    "ids": ids_str,
                }

        if not raw_rows:
            return [], {}

        # Count frequencies per sample
        seq_counter = Counter(sequences)
        final_rows = []
        for r in raw_rows:
            seq = r[3]
            final_rows.append([r[0], r[1], r[2], seq, seq_counter[seq]])

        return final_rows, orfs_meta

    @staticmethod
    def load_cohort_metadata(file_path: Union[str, Path]) -> List[Dict[str, any]]:
        """
        Loads metadata from any TSV or CSV format (runs_metadata.tsv, pilot CSV, or legacy metadata.tsv).
        """
        p = Path(file_path)
        if not p.exists():
            raise FileNotFoundError(f"Metadata file not found: {p}")

        delimiter = "\t" if p.suffix.lower() == ".tsv" else ","
        records = []
        with open(p, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                clean_row = {str(k).strip(): str(v).strip() for k, v in row.items() if k is not None}
                # Find genome_id
                genome_id = (
                    clean_row.get("Genome_ID")
                    or clean_row.get("Genome ID")
                    or clean_row.get("Genome.ID")
                    or clean_row.get("ID")
                    or ""
                ).replace(".tsv", "")

                # Find sra_ids
                sra_str = clean_row.get("SRA Accession") or clean_row.get("Run") or clean_row.get("ID") or ""
                sra_ids = [s.strip().replace(".tsv", "") for s in sra_str.replace(";", ",").split(",") if s.strip()]

                # Find label
                label = (
                    clean_row.get("Label")
                    or clean_row.get("CIP")
                    or clean_row.get("Real_CIP")
                    or clean_row.get("Ciprofloxacin")
                    or clean_row.get("MEM")
                    or ""
                ).upper()

                if genome_id:
                    records.append({
                        "genome_id": genome_id,
                        "sra_ids": sra_ids,
                        "label": label,
                    })
        return records

    def build_dataset(self) -> Tuple[int, int]:
        """
        Main execution method:
        - Reads metadata cohort.
        - Processes all Bakta JSON files.
        - Generates complete_tables and orfs TSVs.
        - Consolidates all_ORFs.json.xz.
        - Exports metadata.tsv.
        """
        logging.info(f"Loading metadata from: {self.metadata_file}")
        records = self.load_cohort_metadata(self.metadata_file)
        logging.info(f"Loaded {len(records)} samples from metadata.")

        all_orfs_catalog: Dict[str, dict] = {}
        processed_metadata_rows: List[Dict[str, str]] = []
        success_count = 0
        missing_count = 0

        for sample in tqdm(records, desc="Processing Bakta annotations"):
            gid = sample["genome_id"]
            sra_ids = sample["sra_ids"]
            label_val = sample["label"]

            json_path = self.find_json_file(gid, sra_ids)
            if not json_path:
                logging.warning(f"[{gid}] Bakta JSON not found in '{self.jsons_dir}'. Skipping.")
                missing_count += 1
                continue

            rows, orfs_meta = self.process_single_json(json_path)
            if not rows:
                logging.warning(f"[{gid}] No valid CDS features extracted from {json_path.name}.")
                missing_count += 1
                continue

            # Merge into global catalog
            for seq, meta in orfs_meta.items():
                if seq not in all_orfs_catalog:
                    all_orfs_catalog[seq] = meta

            # Choose identifier for output filenames
            if self.use_sra_id and sra_ids:
                sample_stem = sra_ids[0]
            else:
                sample_stem = gid

            table_filename = f"{sample_stem}.tsv"

            # 1. Write intermediate complete table [locus, type, product, sequence, count]
            complete_table_path = self.complete_tables_dir / table_filename
            with open(complete_table_path, "w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(["locus", "type", "product", "sequence", "count"])
                writer.writerows(rows)

            # 2. Write DeepRC input table [orf, templates]
            orf_table_path = self.orfs_dir / table_filename
            with open(orf_table_path, "w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(["orf", "templates"])
                for r in rows:
                    writer.writerow([r[3], r[4]])

            # 3. Add to metadata rows
            # Format label as 'R' or 'S'
            if label_val in ["RESISTANT", "R"]:
                canonical_label = "R"
            elif label_val in ["SUSCEPTIBLE", "S"]:
                canonical_label = "S"
            else:
                canonical_label = label_val

            processed_metadata_rows.append({
                "ID": table_filename,
                self.antibiotic_name: canonical_label,
                f"Real_{self.antibiotic_name}": canonical_label,
            })

            success_count += 1

        # 4. Write database/metadata.tsv
        with open(self.metadata_output, "w", encoding="utf-8", newline="") as f:
            fieldnames = ["ID", self.antibiotic_name, f"Real_{self.antibiotic_name}"]
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(processed_metadata_rows)
        logging.info(f"Saved DeepMDC metadata to: {self.metadata_output} ({len(processed_metadata_rows)} samples)")

        # 5. Write database/all_ORFs.json.xz
        logging.info(f"Exporting {len(all_orfs_catalog)} unique ORFs to: {self.all_orfs_output}...")
        if str(self.all_orfs_output).endswith(".xz"):
            with lzma.open(self.all_orfs_output, "wt", encoding="utf-8") as out:
                json.dump(all_orfs_catalog, out)
        else:
            with open(self.all_orfs_output, "w", encoding="utf-8") as out:
                json.dump(all_orfs_catalog, out, indent=2)

        logging.info(f"Successfully generated all_ORFs catalog with {len(all_orfs_catalog)} unique sequences.")
        return success_count, missing_count


def main():
    parser = argparse.ArgumentParser(
        description="DeepMDC ORF Dataset Builder: processes Bakta annotations into orfs/*.tsv and all_ORFs.json.xz.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--metadata", "-m",
        type=str,
        default="database/runs_metadata.tsv",
        help="Path to cohort metadata (TSV or CSV).",
    )
    parser.add_argument(
        "--jsons-dir", "-j",
        type=str,
        default="jsons",
        help="Directory containing Bakta output JSON files (.json or .json.xz).",
    )
    parser.add_argument(
        "--complete-tables-dir",
        type=str,
        default="database/complete_tables",
        help="Directory to save intermediate detailed complete tables.",
    )
    parser.add_argument(
        "--orfs-dir", "-o",
        type=str,
        default="database/orfs",
        help="Directory to save final DeepRC ORF input tables (columns: orf, templates).",
    )
    parser.add_argument(
        "--all-orfs-output", "-a",
        type=str,
        default="database/all_ORFs.json.xz",
        help="Path to save the deduplicated global biological dictionary (.json or .json.xz).",
    )
    parser.add_argument(
        "--metadata-output",
        type=str,
        default="database/metadata.tsv",
        help="Path to save the final DeepMDC metadata.tsv.",
    )
    parser.add_argument(
        "--antibiotic",
        type=str,
        default="CIP",
        help="Name of target antibiotic column (e.g. CIP, MEM).",
    )
    parser.add_argument(
        "--filter-orfs",
        action="store_true",
        help="Filter out hypothetical proteins and sequences outside 40-2000 aa range.",
    )
    parser.add_argument(
        "--use-sra-id",
        action="store_true",
        help="Use SRA run accession (e.g. ERR4035678.tsv) instead of Genome ID (562.100134.tsv) as filename.",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    builder = OrfDatasetBuilder(
        metadata_file=args.metadata,
        jsons_dir=args.jsons_dir,
        complete_tables_dir=args.complete_tables_dir,
        orfs_dir=args.orfs_dir,
        all_orfs_output=args.all_orfs_output,
        metadata_output=args.metadata_output,
        antibiotic_name=args.antibiotic,
        filter_orfs=args.filter_orfs,
        use_sra_id=args.use_sra_id,
    )

    success, missing = builder.build_dataset()

    print("\n" + "=" * 60)
    print("DEEPMDC ORF DATASET GENERATION SUMMARY")
    print("=" * 60)
    print(f"Successfully processed samples: {success}")
    print(f"Missing / skipped samples: {missing}")
    print(f"ORF tables saved to: {args.orfs_dir}/")
    print(f"Complete tables saved to: {args.complete_tables_dir}/")
    print(f"DeepMDC metadata saved to: {args.metadata_output}")
    print(f"Global biological catalog: {args.all_orfs_output}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
