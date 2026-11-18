"""
DeepMDC ORF & Non-ORF Dataset Generator.

Processes Bakta JSON annotation files (.json or .json.xz) into the exact format
expected by DeepMDC:
1. ORFs (Coding sequences - CDS):
   - database/complete_tables/<ID>.tsv: [locus, type, product, sequence, count]
   - database/orfs/<ID>.tsv: [orf, templates]
   - database/all_ORFs.json.xz: Deduplicated amino acid dictionary
2. Non-ORFs (ncRNAs, tRNAs, rRNAs, CRISPRs, regulatory elements):
   - database/complete_tables_non_orfs/<ID>.tsv: [locus, type, product, sequence, count]
   - database/non_orfs/<ID>.tsv: [orf, templates]
   - database/all_non_orfs.json.xz: Deduplicated nucleotide dictionary
3. Metadata:
   - database/metadata.tsv & database/metadata_non_orfs.tsv
"""

import argparse
import logging
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable

_pkg_root = Path(__file__).resolve().parent.parent
if str(_pkg_root) not in sys.path:
    sys.path.insert(0, str(_pkg_root))

from deepmdc_rawdata.io_utils import (
    print_summary_banner,
    read_json_any,
    write_dict_tsv,
    write_json_any,
    write_tsv,
)
from deepmdc_rawdata.metadata import SraMetadata


@dataclass
class OrfExtractionResult:
    """Stores metrics for a single sample's extracted ORFs and non-ORFs."""
    sample_id: str
    total_orfs: int
    unique_orfs: int
    complete_table_path: Path
    orf_table_path: Path
    total_non_orfs: int = 0
    unique_non_orfs: int = 0
    complete_table_non_orf_path: Optional[Path] = None
    non_orf_table_path: Optional[Path] = None


class OrfDatasetBuilder:
    """Builds the DeepMDC ORF and non-ORF datasets and consolidated biological dictionaries."""

    def __init__(
        self,
        metadata_file: Union[str, Path],
        jsons_dir: Union[str, Path] = "jsons",
        complete_tables_dir: Union[str, Path] = "database/complete_tables",
        orfs_dir: Union[str, Path] = "database/orfs",
        all_orfs_output: Union[str, Path] = "database/all_ORFs.json.xz",
        complete_tables_non_orfs_dir: Union[str, Path] = "database/complete_tables_non_orfs",
        non_orfs_dir: Union[str, Path] = "database/non_orfs",
        all_non_orfs_output: Union[str, Path] = "database/all_non_orfs.json.xz",
        metadata_output: Union[str, Path] = "database/metadata.tsv",
        metadata_non_orfs_output: Optional[Union[str, Path]] = "database/metadata_non_orfs.tsv",
        extract_non_orfs: bool = True,
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
        self.complete_tables_non_orfs_dir = Path(complete_tables_non_orfs_dir)
        self.non_orfs_dir = Path(non_orfs_dir)
        self.all_non_orfs_output = Path(all_non_orfs_output)
        self.metadata_output = Path(metadata_output)
        self.metadata_non_orfs_output = Path(metadata_non_orfs_output) if metadata_non_orfs_output else None
        self.extract_non_orfs = extract_non_orfs
        self.antibiotic_name = antibiotic_name.upper()
        self.filter_orfs = filter_orfs
        self.use_sra_id = use_sra_id
        self.min_aa_len = min_aa_len
        self.max_aa_len = max_aa_len

        self.complete_tables_dir.mkdir(parents=True, exist_ok=True)
        self.orfs_dir.mkdir(parents=True, exist_ok=True)
        if self.extract_non_orfs:
            self.complete_tables_non_orfs_dir.mkdir(parents=True, exist_ok=True)
            self.non_orfs_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def read_compressed_json(json_path: Union[str, Path]) -> dict:
        """Reads a plain (.json) or LZMA-compressed (.json.xz) JSON file."""
        return read_json_any(json_path)

    def find_json_file(self, genome_id: str, sra_ids: List[str]) -> Optional[Path]:
        """Locates the Bakta JSON file matching either genome_id or any of its SRA run IDs."""
        for stem in [genome_id, *sra_ids]:
            for ext in (".json.xz", ".json"):
                cand = self.jsons_dir / f"{stem}{ext}"
                if cand.exists() and cand.is_file():
                    return cand
        return None

    def process_single_json(
        self, json_file: Path
    ) -> Tuple[
        List[List[Union[str, int]]],
        Dict[str, dict],
        List[List[Union[str, int]]],
        Dict[str, dict],
    ]:
        """Extracts both CDS amino acid sequences (ORFs) and non-CDS features (non-ORFs)."""
        data = read_json_any(json_file)
        features = data.get("features", [])

        orf_seqs, raw_orf_rows, orfs_meta = [], [], {}
        non_orf_seqs, raw_non_orf_rows, non_orfs_meta = [], [], {}

        for feat in features:
            raw_type = str(feat.get("type", "")).strip()
            feat_type = raw_type.lower()
            product = str(feat.get("product", "")).strip()
            locus = str(feat.get("locus", feat.get("id", ""))).strip()
            gene = str(feat["gene"]).strip() if feat.get("gene") else None
            db_xrefs = feat.get("db_xrefs", [])
            ids_str = ", ".join(str(x) for x in db_xrefs) if isinstance(db_xrefs, list) else str(db_xrefs or "")

            if feat_type == "cds":
                seq = str(feat.get("aa", "")).strip().upper()
                if not seq:
                    continue
                if self.filter_orfs:
                    if product.lower() == "hypothetical protein" or not (self.min_aa_len <= len(seq) <= self.max_aa_len):
                        continue
                orf_seqs.append(seq)
                raw_orf_rows.append([locus, "cds", product, seq])
                if seq not in orfs_meta:
                    orfs_meta[seq] = {"gene": gene, "product": product, "ids": ids_str}
            else:
                nt_seq = str(feat.get("nt", "")).strip().upper()
                if not nt_seq:
                    continue
                non_orf_seqs.append(nt_seq)
                raw_non_orf_rows.append([locus, raw_type, product, nt_seq])
                if nt_seq not in non_orfs_meta:
                    non_orfs_meta[nt_seq] = {"type": raw_type, "gene": gene, "product": product, "ids": ids_str}

        orf_counts = Counter(orf_seqs)
        final_orf_rows = [[r[0], r[1], r[2], r[3], orf_counts[r[3]]] for r in raw_orf_rows]

        non_orf_counts = Counter(non_orf_seqs)
        final_non_orf_rows = [[r[0], r[1], r[2], r[3], non_orf_counts[r[3]]] for r in raw_non_orf_rows]

        return final_orf_rows, orfs_meta, final_non_orf_rows, non_orfs_meta

    def build_dataset(self) -> Tuple[int, int]:
        """Processes all Bakta JSON files and writes ORF/non-ORF tables, catalogs, and metadata."""
        metadata = SraMetadata.load(self.metadata_file)
        logging.info(f"Loaded {len(metadata.samples)} samples from {self.metadata_file}")

        all_orfs_catalog: Dict[str, dict] = {}
        all_non_orfs_catalog: Dict[str, dict] = {}
        processed_metadata_rows: List[Dict[str, str]] = []
        success_count, missing_count = 0, 0

        for sample in tqdm(metadata.samples, desc="Processing Bakta annotations"):
            gid = sample.genome_id
            sra_ids = sample.sra_accessions
            label_val = sample.label.upper()

            json_path = self.find_json_file(gid, sra_ids)
            if not json_path:
                logging.warning(f"[{gid}] Bakta JSON not found in '{self.jsons_dir}'. Skipping.")
                missing_count += 1
                continue

            orf_rows, orfs_meta, non_orf_rows, non_orfs_meta = self.process_single_json(json_path)
            if not orf_rows and not non_orf_rows:
                logging.warning(f"[{gid}] No valid features extracted from {json_path.name}.")
                missing_count += 1
                continue

            for seq, meta in orfs_meta.items():
                all_orfs_catalog.setdefault(seq, meta)
            if self.extract_non_orfs:
                for seq, meta in non_orfs_meta.items():
                    all_non_orfs_catalog.setdefault(seq, meta)

            sample_stem = sra_ids[0] if (self.use_sra_id and sra_ids) else gid
            table_filename = f"{sample_stem}.tsv"

            if orf_rows:
                write_tsv(
                    self.complete_tables_dir / table_filename,
                    ["locus", "type", "product", "sequence", "count"],
                    orf_rows,
                )
                write_tsv(
                    self.orfs_dir / table_filename,
                    ["orf", "templates"],
                    [[r[3], r[4]] for r in orf_rows],
                )

            if self.extract_non_orfs and non_orf_rows:
                write_tsv(
                    self.complete_tables_non_orfs_dir / table_filename,
                    ["locus", "type", "product", "sequence", "count"],
                    non_orf_rows,
                )
                write_tsv(
                    self.non_orfs_dir / table_filename,
                    ["orf", "templates"],
                    [[r[3], r[4]] for r in non_orf_rows],
                )

            canonical_label = "R" if label_val in ("RESISTANT", "R") else ("S" if label_val in ("SUSCEPTIBLE", "S") else label_val)
            processed_metadata_rows.append({
                "ID": table_filename,
                self.antibiotic_name: canonical_label,
                f"Real_{self.antibiotic_name}": canonical_label,
            })
            success_count += 1

        fieldnames = ["ID", self.antibiotic_name, f"Real_{self.antibiotic_name}"]
        write_dict_tsv(self.metadata_output, fieldnames, processed_metadata_rows)
        if self.extract_non_orfs and self.metadata_non_orfs_output:
            write_dict_tsv(self.metadata_non_orfs_output, fieldnames, processed_metadata_rows)

        write_json_any(self.all_orfs_output, all_orfs_catalog)
        if self.extract_non_orfs:
            write_json_any(self.all_non_orfs_output, all_non_orfs_catalog)

        return success_count, missing_count


def main():
    parser = argparse.ArgumentParser(
        description="DeepMDC ORF & Non-ORF Dataset Builder.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metadata", "-m", type=str, default="database/runs_metadata.tsv", help="Path to cohort metadata.")
    parser.add_argument("--jsons-dir", "-j", type=str, default="jsons", help="Directory containing Bakta JSON files.")
    parser.add_argument("--complete-tables-dir", type=str, default="database/complete_tables", help="ORF complete tables dir.")
    parser.add_argument("--orfs-dir", "-o", type=str, default="database/orfs", help="Final DeepRC ORF input tables dir.")
    parser.add_argument("--all-orfs-output", "-a", type=str, default="database/all_ORFs.json.xz", help="Global ORF catalog path.")
    parser.add_argument("--complete-tables-non-orfs-dir", type=str, default="database/complete_tables_non_orfs", help="Non-ORF complete tables dir.")
    parser.add_argument("--non-orfs-dir", type=str, default="database/non_orfs", help="Final DeepRC non-ORF tables dir.")
    parser.add_argument("--all-non-orfs-output", type=str, default="database/all_non_orfs.json.xz", help="Global non-ORF catalog path.")
    parser.add_argument("--skip-non-orfs", action="store_true", help="Skip extraction of non-ORF features.")
    parser.add_argument("--metadata-output", type=str, default="database/metadata.tsv", help="Output metadata.tsv path.")
    parser.add_argument("--metadata-non-orfs-output", type=str, default="database/metadata_non_orfs.tsv", help="Output non-ORF metadata.tsv path.")
    parser.add_argument("--antibiotic", type=str, default="CIP", help="Target antibiotic column name.")
    parser.add_argument("--filter-orfs", action="store_true", help="Filter out hypothetical proteins and sequences outside 40-2000 aa.")
    parser.add_argument("--use-sra-id", action="store_true", help="Use SRA run accession instead of Genome ID as filename.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    builder = OrfDatasetBuilder(
        metadata_file=args.metadata,
        jsons_dir=args.jsons_dir,
        complete_tables_dir=args.complete_tables_dir,
        orfs_dir=args.orfs_dir,
        all_orfs_output=args.all_orfs_output,
        complete_tables_non_orfs_dir=args.complete_tables_non_orfs_dir,
        non_orfs_dir=args.non_orfs_dir,
        all_non_orfs_output=args.all_non_orfs_output,
        metadata_output=args.metadata_output,
        metadata_non_orfs_output=args.metadata_non_orfs_output,
        extract_non_orfs=not args.skip_non_orfs,
        antibiotic_name=args.antibiotic,
        filter_orfs=args.filter_orfs,
        use_sra_id=args.use_sra_id,
    )

    success, missing = builder.build_dataset()

    summary_items = [
        ("Successfully processed samples", success),
        ("Missing / skipped samples", missing),
        ("ORF tables saved to", f"{args.orfs_dir}/"),
        ("ORF complete tables saved to", f"{args.complete_tables_dir}/"),
        ("Global ORF catalog", args.all_orfs_output),
    ]
    if not args.skip_non_orfs:
        summary_items.extend([
            ("Non-ORF tables saved to", f"{args.non_orfs_dir}/"),
            ("Non-ORF complete tables saved to", f"{args.complete_tables_non_orfs_dir}/"),
            ("Global non-ORF catalog", args.all_non_orfs_output),
        ])
    summary_items.append(("DeepMDC metadata saved to", args.metadata_output))
    print_summary_banner("DEEPMDC ORF & NON-ORF DATASET GENERATION SUMMARY", summary_items)


if __name__ == "__main__":
    main()
