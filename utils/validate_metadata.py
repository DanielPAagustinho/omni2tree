#!/usr/bin/env python3
"""
Validate Omni2tree Step 3 metadata before running expensive steps.

Checks:
1) Metadata CSV structure (header + type row + data rows)
2) Required columns (label, accession)
3) Allowed type declarations
4) label/accession constraints (non-empty; accession one-per-row; unique label/accession)
5) Label collision risk after output sanitization
6) Every reference label in five_letter_taxon.tsv is present in metadata
   after alphanumeric cleanup (same criterion used in step1/step3 matching)
7) Every readset with consensus (O2T_RESULTS/*_all_cov.txt) is present in metadata label column
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple


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


VALID_TYPES = {"character", "date", "numeric", "integer"}


@dataclass
class MetadataInput:
    header: List[str]
    types: List[str]
    rows: List[Dict[str, str]]
    label_col: str
    accession_col: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Validate step3 metadata and reference coverage against five_letter_taxon.tsv"
    )
    p.add_argument("-m", "--metadata", required=True, type=Path, help="Input metadata CSV")
    p.add_argument(
        "--five_letter",
        required=True,
        type=Path,
        help="five_letter_taxon.tsv from step1",
    )
    p.add_argument(
        "--o2t_results",
        default=Path("O2T_RESULTS"),
        type=Path,
        help="Read2Tree output directory containing *_all_cov.txt files (default: ./O2T_RESULTS)",
    )
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

    bad_types = [t for t in types if t not in VALID_TYPES]
    if bad_types:
        raise ValueError(
            f"Invalid metadata column types in second row: {', '.join(sorted(set(bad_types)))} "
            f"(allowed: {', '.join(sorted(VALID_TYPES))})"
        )

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

        label = row[label_col]
        if label == "":
            raise ValueError(f"Metadata row {i} has empty label")
        if "," in label:
            raise ValueError(
                f"Metadata row {i} label contains a comma; labels must not contain commas"
            )

        accession = row[accession_col]
        if accession == "":
            raise ValueError(f"Metadata row {i} has empty accession")
        if "," in accession:
            raise ValueError(
                f"Metadata row {i} accession contains a comma; exactly one accession is required per row"
            )

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
    seen_codes: set[str] = set()

    with path.open(encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Invalid five_letter line {line_no}: expected 2 tab-separated columns")

            taxon, code = parts[0].strip(), parts[1].strip()
            if code in seen_codes:
                raise ValueError(f"Duplicated five-letter code in {path}: {code}")
            seen_codes.add(code)

            label_key = clean_alnum(taxon)
            if label_key in labelkey_to_code:
                raise ValueError(
                    "Colliding taxon names in five_letter_taxon.tsv after alphanumeric cleanup: "
                    f"'{taxon}' and '{code_to_taxon[labelkey_to_code[label_key]]}'"
                )

            labelkey_to_code[label_key] = code
            code_to_taxon[code] = taxon

    log_info(f"Loaded {len(code_to_taxon)} entries from {path}")
    return labelkey_to_code, code_to_taxon


def validate_reference_coverage(meta: MetadataInput, labelkey_to_code: Dict[str, str], code_to_taxon: Dict[str, str]) -> None:
    metadata_label_keys: set[str] = set()
    for row in meta.rows:
        metadata_label_keys.add(clean_alnum(row[meta.label_col]))

    missing_keys = sorted(set(labelkey_to_code) - metadata_label_keys)
    if missing_keys:
        missing_refs = [code_to_taxon[labelkey_to_code[k]] for k in missing_keys]
        raise ValueError(
            "Metadata is missing reference label(s) required by five_letter_taxon.tsv "
            "(comparison after alphanumeric cleanup): " + ", ".join(missing_refs)
        )

    log_info("Validated metadata coverage for all references in five_letter_taxon.tsv")


def load_readsets_from_cov(o2t_results: Path) -> List[str]:
    if not o2t_results.exists():
        raise FileNotFoundError(f"O2T_RESULTS directory not found: {o2t_results}")
    if not o2t_results.is_dir():
        raise ValueError(f"O2T_RESULTS path is not a directory: {o2t_results}")

    suffix = "_all_cov.txt"
    names = sorted(
        {
            p.name[:-len(suffix)]
            for p in o2t_results.glob(f"*{suffix}")
            if p.is_file() and p.name.endswith(suffix) and p.name != suffix
        }
    )
    if not names:
        raise ValueError(f"No '*{suffix}' files found in {o2t_results}")

    log_info(f"Detected {len(names)} consensus readset(s) from {o2t_results}")
    return names


def validate_readset_presence(meta: MetadataInput, readset_names: List[str]) -> None:
    labels = {row[meta.label_col] for row in meta.rows}
    missing = [name for name in readset_names if name not in labels]
    if missing:
        raise ValueError(
            "Metadata label column is missing readset name(s) present in O2T_RESULTS/*_all_cov.txt: "
            + ", ".join(missing)
        )
    log_info("Validated metadata label coverage for all readsets with consensus (_all_cov.txt)")


def validate_output_constraints(meta: MetadataInput, labelkey_to_code: Dict[str, str]) -> None:
    seen_sample_ids: set[str] = set()
    seen_raw_labels: set[str] = set()
    seen_label_keys: Dict[str, str] = {}
    ref_count = 0
    read_count = 0

    for idx, row in enumerate(meta.rows, start=1):
        accession = row[meta.accession_col]
        raw_label = row[meta.label_col]

        if raw_label in seen_raw_labels:
            raise ValueError(f"Duplicated label in metadata: {raw_label}")
        seen_raw_labels.add(raw_label)

        if accession in seen_sample_ids:
            raise ValueError(f"Duplicated accession/sample_id in metadata: {accession}")
        seen_sample_ids.add(accession)

        label_key = clean_alnum(raw_label)
        label_safe = tree_label_sanitize(raw_label)
        label_final = re.sub(r"[^A-Za-z0-9_]", "", label_safe) or "NA"
        collision_key = clean_alnum(label_final)
        if collision_key in seen_label_keys:
            raise ValueError(
                "Label collision after alphanumeric normalization: "
                f"'{raw_label}' collides with '{seen_label_keys[collision_key]}'"
            )
        seen_label_keys[collision_key] = raw_label

        if label_key in labelkey_to_code:
            ref_count += 1
        else:
            read_count += 1

    log_info(
        "Validated metadata output constraints: "
        f"{len(meta.rows)} rows ({ref_count} references, {read_count} readsets)"
    )


def main() -> None:
    args = parse_args()

    try:
        meta = parse_metadata_csv(args.metadata)
        labelkey_to_code, code_to_taxon = load_five_letter(args.five_letter)
        readset_names = load_readsets_from_cov(args.o2t_results)
        validate_reference_coverage(meta, labelkey_to_code, code_to_taxon)
        validate_readset_presence(meta, readset_names)
        validate_output_constraints(meta, labelkey_to_code)
    except Exception as exc:
        log_error(str(exc))
        sys.exit(1)

    log_info("validate_metadata.py completed successfully")


if __name__ == "__main__":
    main()
