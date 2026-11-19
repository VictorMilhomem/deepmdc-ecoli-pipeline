"""
Cohort metadata models and parsers for Escherichia coli AMR datasets.
"""

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

from .io_utils import get_project_root, resolve_path, write_tsv

PILOT_METADATA_REL = "database/ecoli_ciprofloxacin_RS_pilot_50R_50S.csv"
CANDIDATES_METADATA_REL = "database/ecoli_ciprofloxacin_RS_candidates_with_metadata.csv"


def resolve_metadata_path(
    file_arg: Optional[str] = None,
    use_all: bool = False,
    default_rel: str = PILOT_METADATA_REL,
) -> Path:
    """Resolves the target metadata path based on CLI flags (--file, --all, or default pilot)."""
    if use_all:
        return get_project_root() / CANDIDATES_METADATA_REL
    if file_arg:
        p = Path(file_arg)
        return p if p.exists() or p.is_absolute() else resolve_path(file_arg)
    return resolve_path(default_rel)


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
    ids: List[str] = field(default_factory=list)
    samples: List[SampleRecord] = field(default_factory=list)
    run_to_sample: Dict[str, SampleRecord] = field(default_factory=dict)
    mems: List[str] = field(default_factory=list)
    real_mems: List[str] = field(default_factory=list)

    @classmethod
    def from_csv(cls, csv_file: Union[str, Path], delimiter: str = ";") -> "SraMetadata":
        """
        Parses an E. coli metadata CSV file (pilot or full candidates CSV).
        Handles UTF-8 with BOM (utf-8-sig), auto-detects ';' vs ',', treats Genome.ID as string,
        and expands multiple comma-separated SRA runs per isolate.
        """
        metadata = cls()
        csv_path = Path(csv_file)
        if not csv_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {csv_path}")

        with open(csv_path, mode="r", encoding="utf-8-sig") as f:
            first_line = f.readline()
            f.seek(0)
            actual_delim = "," if (delimiter == ";" and ";" not in first_line and "," in first_line) else delimiter

            for row in csv.DictReader(f, delimiter=actual_delim):
                clean = {k.strip(): (v.strip() if v is not None else "") for k, v in row.items() if k is not None}

                genome_name = clean.get("Genome.Name") or clean.get("Genome Name") or clean.get("Genome_Name") or ""
                sra_raw = clean.get("SRA.Accession") or clean.get("SRA Accession") or clean.get("Run") or clean.get("ID") or ""
                cipro = clean.get("ciprofloxacin") or clean.get("Ciprofloxacin") or ""
                label = clean.get("label") or clean.get("Label") or clean.get("CIP") or clean.get("Real_CIP") or ""
                genome_id = str(
                    clean.get("Genome.ID") or clean.get("Genome ID") or clean.get("Genome_ID") or clean.get("ID") or ""
                ).replace(".tsv", "").strip()

                runs = [r.strip().replace(".tsv", "") for r in sra_raw.replace(";", ",").split(",") if r.strip()]

                sample = SampleRecord(
                    genome_name=genome_name,
                    sra_accessions=runs,
                    ciprofloxacin=cipro,
                    label=label,
                    genome_id=genome_id,
                    biosample_accession=clean.get("BioSample.Accession") or clean.get("BioSample Accession") or clean.get("BioSample_Accession") or "",
                    assembly_accession=clean.get("Assembly.Accession") or clean.get("Assembly Accession") or clean.get("Assembly_Accession") or "",
                    genome_quality=clean.get("Genome.Quality") or clean.get("Genome Quality") or "",
                    contigs=clean.get("Contigs", ""),
                    contig_l50=clean.get("Contig.L50") or clean.get("Contig L50") or "",
                    contig_n50=clean.get("Contig.N50") or clean.get("Contig N50") or "",
                    checkm_completeness=clean.get("CheckM.Completeness") or clean.get("CheckM Completeness") or "",
                    checkm_contamination=clean.get("CheckM.Contamination") or clean.get("CheckM Contamination") or "",
                    mlst=clean.get("MLST", ""),
                    host_name=clean.get("Host.Name") or clean.get("Host Name") or "",
                    isolation_source=clean.get("Isolation.Source") or clean.get("Isolation Source") or "",
                    geographic_location=clean.get("Geographic.Location") or clean.get("Geographic Location") or "",
                    extra_fields=clean,
                )
                metadata.samples.append(sample)

                for run_id in runs:
                    metadata.ids.append(run_id)
                    metadata.run_to_sample[run_id] = sample
                    metadata.mems.append(label)
                    metadata.real_mems.append(label)

        return metadata

    @classmethod
    def from_tsv(cls, tsv_file: Union[str, Path]) -> "SraMetadata":
        """Parses a TSV metadata file (runs_metadata.tsv or legacy metadata.tsv)."""
        metadata = cls()
        tsv_path = Path(tsv_file)
        if not tsv_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {tsv_path}")

        with open(tsv_path, mode="r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if not header:
                return metadata

            clean_header = [h.strip() for h in header]
            if any(k in clean_header for k in ("Genome_ID", "Genome.ID", "Genome ID")):
                f.seek(0)
                seen_genomes: Dict[str, SampleRecord] = {}
                for row in csv.DictReader(f, delimiter="\t"):
                    clean = {str(k).strip(): str(v).strip() for k, v in row.items() if k is not None}
                    run_id = (clean.get("ID") or clean.get("Run") or "").replace(".tsv", "")
                    genome_id = (
                        clean.get("Genome_ID") or clean.get("Genome.ID") or clean.get("Genome ID") or run_id
                    ).replace(".tsv", "")
                    label = clean.get("Label") or clean.get("CIP") or clean.get("Real_CIP") or clean.get("Ciprofloxacin") or clean.get("MEM") or ""

                    if not genome_id and not run_id:
                        continue

                    if genome_id in seen_genomes:
                        sample = seen_genomes[genome_id]
                        if run_id and run_id not in sample.sra_accessions:
                            sample.sra_accessions.append(run_id)
                    else:
                        sample = SampleRecord(
                            genome_name=clean.get("Genome_Name") or clean.get("Genome.Name") or "",
                            sra_accessions=[run_id] if run_id else [],
                            ciprofloxacin=clean.get("Ciprofloxacin") or clean.get("ciprofloxacin") or "",
                            label=label,
                            genome_id=genome_id,
                            biosample_accession=clean.get("BioSample_Accession") or "",
                            assembly_accession=clean.get("Assembly_Accession") or "",
                        )
                        seen_genomes[genome_id] = sample
                        metadata.samples.append(sample)

                    if run_id:
                        metadata.ids.append(run_id)
                        metadata.run_to_sample[run_id] = sample
                        metadata.mems.append(label)
                        metadata.real_mems.append(label)
            else:
                for row in reader:
                    if len(row) >= 2:
                        srr_id = row[0].replace(".tsv", "").strip()
                        label = row[1].strip()
                        real_label = row[2].strip() if len(row) >= 3 else label
                        sample = SampleRecord(
                            genome_name="",
                            sra_accessions=[srr_id],
                            ciprofloxacin=label,
                            label=label,
                            genome_id=srr_id,
                        )
                        metadata.samples.append(sample)
                        metadata.ids.append(srr_id)
                        metadata.run_to_sample[srr_id] = sample
                        metadata.mems.append(label)
                        metadata.real_mems.append(real_label)
        return metadata

    @classmethod
    def load(cls, file_path: Union[str, Path]) -> "SraMetadata":
        """Auto-detects file type (.csv or .tsv) and loads metadata."""
        path = Path(file_path)
        return cls.from_tsv(path) if path.suffix.lower() == ".tsv" else cls.from_csv(path)

    def summary(self) -> Dict[str, int]:
        """Returns summary statistics for the loaded metadata."""
        return {
            "total_isolates": len(self.samples),
            "total_runs": len(self.ids),
            "resistant": sum(1 for s in self.samples if s.label.upper() in ("R", "RESISTANT")),
            "susceptible": sum(1 for s in self.samples if s.label.upper() in ("S", "SUSCEPTIBLE")),
            "multi_run_isolates": sum(1 for s in self.samples if len(s.sra_accessions) > 1),
        }

    def export_run_manifest_tsv(self, output_tsv: Union[str, Path]) -> None:
        """Exports a TSV manifest with one row per SRA run."""
        header = [
            "ID",
            "Genome_ID",
            "Genome_Name",
            "Ciprofloxacin",
            "Label",
            "BioSample_Accession",
            "Assembly_Accession",
        ]
        rows = []
        for run_id in self.ids:
            s = self.run_to_sample.get(run_id)
            if s:
                rows.append([
                    run_id,
                    s.genome_id,
                    s.genome_name,
                    s.ciprofloxacin,
                    s.label,
                    s.biosample_accession,
                    s.assembly_accession,
                ])
            else:
                rows.append([run_id, "", "", "", "", "", ""])
        write_tsv(output_tsv, header, rows)
