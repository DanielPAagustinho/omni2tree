#!/usr/bin/env python3
"""
Prepare Omni2tree visualization inputs:
1) Rewrite tree leaf labels replacing step1 five-letter codes with real labels.
2) Build a metadata CSV compatible with omni2treeview.py (header + type row).

Notes:
- Input metadata must include a first data row with column types.
- Required metadata columns: label, accession (case-insensitive).
- Output metadata uses sample_id=accession (as requested for step3 notes).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import Phylo


if sys.stdout.isatty():
    RED = "\033[1;31m"
    GREEN = "\033[1;32m"
    YELLOW = "\033[1;33m"
    NC = "\033[0m"
else:
    RED = GREEN = YELLOW = NC = ""


def log_info(message: str) -> None:
    print(f"{GREEN}[{datetime.now():%Y-%m-%d %H:%M:%S}] [INFO]{NC} {message}")


def log_warn(message: str) -> None:
    print(f"{YELLOW}[{datetime.now():%Y-%m-%d %H:%M:%S}] [WARN]{NC} {message}")


def log_error(message: str) -> None:
    print(f"{RED}[{datetime.now():%Y-%m-%d %H:%M:%S}] [ERROR]{NC} {message}", file=sys.stderr)


@dataclass
class MetadataInput:
    header: List[str]
    types: List[str]
    rows: List[Dict[str, str]]
    label_col: str
    accession_col: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare metadata and relabeled tree for omni2treeview"
    )
    p.add_argument("-m", "--metadata", required=True, type=Path, help="Input metadata CSV")
    p.add_argument("--five_letter", default=Path("five_letter_taxon.tsv"), type=Path,
                   help="five_letter_taxon.tsv from step1 (default: ./five_letter_taxon.tsv)")
    p.add_argument("--in_nwk", required=True, type=Path, help="Input Newick tree")
    p.add_argument("--out_nwk", required=True, type=Path, help="Output relabeled Newick tree")
    p.add_argument("--out_meta", required=True, type=Path, help="Output metadata CSV for omni2treeview")
    return p.parse_args()


def clean_alnum(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", text or "")


def tree_label_sanitize(text: str) -> str:
    text = re.sub(r"_(1|2)$", "", text or "")
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "NA"


def parse_metadata_csv(path: Path) -> MetadataInput:
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        rows = list(reader)

    if len(rows) < 3:
        raise ValueError("Metadata must contain header + type row + at least one data row")

    header = [h.strip() for h in rows[0]]
    types = [t.strip().lower() for t in rows[1]]
    if len(header) != len(types):
        raise ValueError("Metadata header and type row have different column counts")
    lower_map = {col.lower(): col for col in header}
    if "label" not in lower_map:
        raise ValueError("Metadata must contain a 'label' column")
    if "accession" not in lower_map:
        raise ValueError("Metadata must contain an 'accession' column")
    label_col = lower_map["label"]
    accession_col = lower_map["accession"]

    data_rows: List[Dict[str, str]] = []
    for i, raw in enumerate(rows[2:], start=3):
        if len(raw) < len(header):
            raw = raw + [""] * (len(header) - len(raw))
        elif len(raw) > len(header):
            raise ValueError(f"Metadata row {i} has more columns than header")
        row = {header[idx]: raw[idx].strip() for idx in range(len(header))}
        if not any(v != "" for v in row.values()):
            continue

        data_rows.append(row)

    if not data_rows:
        raise ValueError("Metadata has no data rows after filtering empty rows")

    log_info(f"Loaded metadata rows: {len(data_rows)}")
    return MetadataInput(header, types, data_rows, label_col, accession_col)


def load_five_letter(path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"five_letter file not found: {path}")

    labelkey_to_code: Dict[str, str] = {}
    code_to_taxon: Dict[str, str] = {}
    with path.open(encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Invalid five_letter line {line_no}: expected 2 tab-separated columns")
            taxon, code = parts[0].strip(), parts[1].strip()
            label_key = clean_alnum(taxon)
            labelkey_to_code[label_key] = code
            code_to_taxon[code] = taxon

    log_info(f"Loaded {len(code_to_taxon)} entries from {path}")
    return labelkey_to_code, code_to_taxon


def build_output_metadata(
    meta: MetadataInput,
    labelkey_to_code: Dict[str, str],
) -> Tuple[List[str], List[str], List[List[str]], Dict[str, str]]:
    extra_cols = [c for c in meta.header if c.lower() not in {"label", "accession"}]
    extra_types = [meta.types[meta.header.index(c)] for c in extra_cols]

    out_header = ["sample_id", "label", "source"] + extra_cols
    out_types = ["character", "character", "character"] + extra_types

    out_rows: List[List[str]] = []
    code_to_preferred_label: Dict[str, str] = {}
    for row in meta.rows:
        raw_label = row[meta.label_col]
        accession = row[meta.accession_col]
        label_key = clean_alnum(raw_label)
        label_safe = tree_label_sanitize(raw_label)
        label_final = re.sub(r"[^A-Za-z0-9_]", "", label_safe) or "NA"

        is_reference = label_key in labelkey_to_code
        source = "Reference" if is_reference else "Readset"
        if is_reference:
            code = labelkey_to_code[label_key]
            code_to_preferred_label[code] = label_final

        extras = [row.get(c, "") or "NA" for c in extra_cols]
        out_rows.append([accession, label_final, source] + extras)

    n_refs = sum(1 for r in out_rows if r[2] == "Reference")
    n_reads = sum(1 for r in out_rows if r[2] == "Readset")
    log_info(f"Prepared visualization metadata rows: {len(out_rows)} ({n_refs} references, {n_reads} readsets)")
    return out_header, out_types, out_rows, code_to_preferred_label


def ensure_unique_tree_labels(labels: List[str]) -> List[str]:
    counts = Counter()
    out: List[str] = []
    for label in labels:
        counts[label] += 1
        if counts[label] == 1:
            out.append(label)
        else:
            out.append(f"{label}_{counts[label]}")
    return out


def rewrite_tree_labels(
    in_nwk: Path,
    out_nwk: Path,
    code_to_taxon: Dict[str, str],
    code_to_preferred_label: Dict[str, str],
) -> int:
    tree = Phylo.read(str(in_nwk), "newick")
    terminals = list(tree.get_terminals())
    if not terminals:
        raise ValueError(f"No terminal nodes found in tree: {in_nwk}")

    final_names: List[str] = []
    replaced_count = 0
    stripped_pair_suffix = 0

    for clade in terminals:
        original = clade.name or ""
        stripped = re.sub(r"_(1|2)$", "", original)
        if stripped != original:
            stripped_pair_suffix += 1

        new_name = None
        tokens = stripped.split("_")
        for tok in tokens:
            if tok in code_to_taxon:
                preferred = code_to_preferred_label.get(tok, tree_label_sanitize(code_to_taxon[tok]))
                new_name = preferred
                replaced_count += 1
                break

        if new_name is None:
            new_name = tree_label_sanitize(stripped)

        final_names.append(new_name)

    unique_names = ensure_unique_tree_labels(final_names)
    for clade, new_name in zip(terminals, unique_names):
        clade.name = new_name

    out_nwk.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(out_nwk), "newick")
    log_info(f"Wrote relabeled tree: {out_nwk}")
    if stripped_pair_suffix > 0:
        log_info(f"Removed _1/_2 suffixes from {stripped_pair_suffix} terminal name(s)")
    log_info(f"Replaced {replaced_count} terminal name(s) using five-letter code mapping")
    return len(terminals)


def write_output_metadata(out_meta: Path, header: List[str], types: List[str], rows: List[List[str]]) -> None:
    out_meta.parent.mkdir(parents=True, exist_ok=True)
    with out_meta.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(header)
        writer.writerow(types)
        writer.writerows(rows)
    log_info(f"Wrote visualization metadata: {out_meta}")


def main() -> None:
    args = parse_args()

    try:
        meta = parse_metadata_csv(args.metadata)
        labelkey_to_code, code_to_taxon = load_five_letter(args.five_letter)
        out_header, out_types, out_rows, code_to_preferred_label = build_output_metadata(meta, labelkey_to_code)
        n_leaves = rewrite_tree_labels(args.in_nwk, args.out_nwk, code_to_taxon, code_to_preferred_label)
        write_output_metadata(args.out_meta, out_header, out_types, out_rows)
    except Exception as exc:
        log_error(str(exc))
        sys.exit(1)

    log_info("prepare_metadata_o2t_view.py completed successfully")
    log_info(f"Tree leaves processed: {n_leaves}")
    log_info(f"Metadata rows written: {len(out_rows)}")


if __name__ == "__main__":
    main()
