"""
Shared I/O, compression (.json / .json.xz), TSV, and path utilities for DeepMDC-RawData.
"""

import csv
import json
import lzma
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union


def get_project_root() -> Path:
    """Returns the repository root directory (two levels up from src/deepmdc_rawdata)."""
    return Path(__file__).resolve().parent.parent.parent


def resolve_path(path_like: Union[str, Path], base_dir: Optional[Path] = None) -> Path:
    """Resolves a path relative to project root (or base_dir) unless already absolute."""
    p = Path(path_like)
    if p.is_absolute():
        return p
    root = base_dir if base_dir is not None else get_project_root()
    return root / p


def read_json_any(json_path: Union[str, Path]) -> dict:
    """Reads a plain (.json) or LZMA-compressed (.json.xz) JSON file."""
    p = Path(json_path)
    if not p.exists():
        return {}
    if str(p).endswith(".xz"):
        with lzma.open(p, "rt", encoding="utf-8") as f:
            return json.load(f)
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json_any(json_path: Union[str, Path], data: dict, indent: Optional[int] = 2) -> None:
    """Writes a dict to a plain (.json) or LZMA-compressed (.json.xz) JSON file."""
    p = Path(json_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if str(p).endswith(".xz"):
        with lzma.open(p, "wt", encoding="utf-8") as out:
            json.dump(data, out)
    else:
        with open(p, "w", encoding="utf-8") as out:
            json.dump(data, out, indent=indent)


def write_tsv(
    output_path: Union[str, Path],
    header: Sequence[str],
    rows: Iterable[Sequence[Any]],
) -> None:
    """Writes a tab-separated file with header, creating parent directories if needed."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def write_dict_tsv(
    output_path: Union[str, Path],
    fieldnames: Sequence[str],
    rows: Iterable[Dict[str, Any]],
) -> None:
    """Writes a list of dicts to a tab-separated file with the given fieldnames."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def print_summary_banner(title: str, items: Sequence[Tuple[str, Any]]) -> None:
    """Prints a formatted CLI summary banner."""
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)
    for label, value in items:
        if value is None or value == "":
            print(label)
        else:
            print(f"{label}: {value}")
    print("=" * 60 + "\n")
