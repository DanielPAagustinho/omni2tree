#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from extract_local_cds_from_gff import (
    ManifestRow,
    clean_alnum,
    make_seqrecords,
    parse_feature_types,
    parse_group_by,
    resolve_path,
)


if sys.stdout.isatty():
    RED = "\033[1;31m"
    GREEN = "\033[1;32m"
    YELLOW = "\033[1;33m"
    NC = "\033[0m"
else:
    RED = GREEN = YELLOW = NC = ""


def log_info(message: str) -> None:
    print(f"{GREEN}[{datetime.now():%Y-%m-%d %H:%M:%S}] [INFO]{NC} {message}")


def log_error(message: str) -> None:
    print(f"{RED}[{datetime.now():%Y-%m-%d %H:%M:%S}] [ERROR]{NC} {message}", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Validate Omni2Tree local assemblies manifest and selected GTF/GFF3 annotations."
    )
    p.add_argument("--manifest", required=True, type=Path, help="CSV with local assemblies")
    p.add_argument("--local-clean-output", required=True, type=Path, help="Filtered local manifest output")
    p.add_argument("--local-taxa-output", required=True, type=Path, help="Cleaned local taxa output")
    p.add_argument("--status-output", required=True, type=Path, help="Validation status output")
    p.add_argument(
        "--feature-types",
        default="CDS",
        help="Comma-separated GTF/GFF3 feature type(s) to extract from column 3 [default: CDS]",
    )
    p.add_argument(
        "--group-by",
        default=None,
        help="Single GTF/GFF3 attribute used to group selected feature rows into one sequence unit",
    )
    p.add_argument(
        "--code-mode",
        required=True,
        choices=("auto", "with-code", "no-code"),
        help="Expected local manifest code mode",
    )
    p.add_argument("--main-clean-file", type=Path, default=None, help="Cleaned NCBI accession CSV")
    return p.parse_args()


def read_manifest_lines(path: Path) -> List[Tuple[int, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"The local assemblies manifest '{path}' does not exist or is empty.")

    lines = []
    with path.open(encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            if line.strip() and not line.lstrip().startswith("#"):
                lines.append((line_no, line.rstrip("\n")))

    if not lines:
        raise ValueError(
            f"The local assemblies manifest '{path}' has no data after removing comments and empty lines."
        )
    return lines


def header_columns(header_line: str) -> List[str]:
    try:
        return [col.strip() for col in next(csv.reader([header_line]))]
    except csv.Error as exc:
        raise ValueError(f"Invalid local assemblies manifest header: {exc}") from exc


def canonical_header(name: str) -> str:
    return clean_alnum(name).lower()


def detect_has_code_column(lines: List[Tuple[int, str]], code_mode: str) -> bool:
    header = header_columns(lines[0][1])
    if code_mode == "with-code":
        return True
    if code_mode == "no-code":
        return False
    if len(header) == 4:
        return True
    if len(header) == 3:
        return False
    raise ValueError(
        "Local assemblies manifest must have 3 columns (taxon,dna,gff) or 4 columns "
        "(taxon,code,dna,gff) when -i/--input is omitted."
    )


def validate_header(header: List[str], has_code_column: bool) -> None:
    normalized = [canonical_header(col) for col in header]
    if has_code_column:
        if len(normalized) != 4:
            raise ValueError("Local assemblies manifest with codes must have 4 columns: taxon,code,dna,gff")
        taxon_col, code_col, dna_col, ann_col = normalized
        if taxon_col not in {"taxon", "taxa", "label", "labels", "strain", "strains", "species"}:
            raise ValueError(f"Invalid first local manifest column: {header[0]}")
        if code_col not in {"code", "codes"}:
            raise ValueError(f"Invalid second local manifest column; expected code(s), found: {header[1]}")
    else:
        if len(normalized) != 3:
            raise ValueError("Local assemblies manifest without codes must have 3 columns: taxon,dna,gff")
        taxon_col, dna_col, ann_col = normalized
        if taxon_col not in {"taxon", "taxa", "label", "labels", "strain", "strains", "species"}:
            raise ValueError(f"Invalid first local manifest column: {header[0]}")

    if dna_col not in {"dna", "fasta", "fa", "fna", "genome", "assembly", "assemblyfasta"}:
        raise ValueError(f"Invalid DNA/FASTA local manifest column: {header[-2]}")
    if ann_col not in {"gff", "gff3", "gtf", "annotation", "annotations"}:
        raise ValueError(f"Invalid annotation local manifest column: {header[-1]}")


def load_manifest(path: Path, has_code_column: bool) -> List[ManifestRow]:
    base_dir = path.resolve().parent
    rows: List[ManifestRow] = []
    seen_taxa: Dict[str, int] = {}
    seen_codes: Dict[str, int] = {}

    with path.open(newline="", encoding="utf-8") as fh:
        filtered = (line for line in fh if line.strip() and not line.lstrip().startswith("#"))
        reader = csv.reader(filtered)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Local assemblies manifest is empty: {path}") from exc

        header = [col.strip() for col in header]
        validate_header(header, has_code_column)

        for row_idx, row in enumerate(reader, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            expected_cols = 4 if has_code_column else 3
            if len(row) != expected_cols:
                raise ValueError(
                    f"Local manifest line {row_idx} has {len(row)} columns; expected {expected_cols}"
                )

            taxon_raw = row[0].strip()
            taxon = clean_alnum(taxon_raw)
            if not taxon:
                raise ValueError(f"Local manifest line {row_idx} has an empty taxon after cleanup")
            if taxon in seen_taxa:
                raise ValueError(
                    f"Duplicated local taxon after alphanumeric cleanup on line {row_idx}: "
                    f"{taxon_raw} (first seen on line {seen_taxa[taxon]})"
                )
            seen_taxa[taxon] = row_idx

            if has_code_column:
                code = row[1].strip()
                dna_raw = row[2]
                ann_raw = row[3]
                if not re.fullmatch(r"[A-Za-z0-9]{5}", code or ""):
                    raise ValueError(
                        f"Invalid local code on line {row_idx}: '{code}'. "
                        "Codes must have exactly 5 alphanumeric characters."
                    )
                if code in seen_codes:
                    raise ValueError(
                        f"Duplicated local code on line {row_idx}: {code} "
                        f"(first seen on line {seen_codes[code]})"
                    )
                seen_codes[code] = row_idx
            else:
                code = None
                dna_raw = row[1]
                ann_raw = row[2]

            dna_path = resolve_path(dna_raw, base_dir)
            annotation_path = resolve_path(ann_raw, base_dir)
            if not dna_path.is_file() or dna_path.stat().st_size == 0:
                raise ValueError(f"DNA FASTA not found or empty for local taxon {taxon_raw}: {dna_path}")
            if not annotation_path.is_file() or annotation_path.stat().st_size == 0:
                raise ValueError(
                    f"GTF/GFF3 annotation not found or empty for local taxon {taxon_raw}: {annotation_path}"
                )

            rows.append(
                ManifestRow(
                    line_no=row_idx,
                    taxon_raw=taxon_raw,
                    taxon=taxon,
                    code=code,
                    dna_path=dna_path,
                    annotation_path=annotation_path,
                )
            )

    if not rows:
        raise ValueError(f"Local assemblies manifest has no data rows: {path}")
    return rows


def write_clean_manifest(rows: List[ManifestRow], has_code_column: bool, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        if has_code_column:
            writer.writerow(["taxon", "code", "dna", "gff"])
            for row in rows:
                writer.writerow([row.taxon_raw, row.code, row.dna_path, row.annotation_path])
        else:
            writer.writerow(["taxon", "dna", "gff"])
            for row in rows:
                writer.writerow([row.taxon_raw, row.dna_path, row.annotation_path])



def read_main_keys(path: Path, has_code_column: bool) -> Tuple[Dict[str, Tuple[int, str]], Dict[str, Tuple[int, str]]]:
    taxa: Dict[str, Tuple[int, str]] = {}
    codes: Dict[str, Tuple[int, str]] = {}

    if path is None:
        return taxa, codes
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Cleaned accession file not found or empty: {path}")

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(line for line in fh if line.strip() and not line.lstrip().startswith("#"))
        try:
            next(reader)
        except StopIteration:
            return taxa, codes

        for row_idx, row in enumerate(reader, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            taxon_raw = row[0].strip()
            taxon = clean_alnum(taxon_raw)
            if taxon:
                if taxon in taxa:
                    first_line, first_raw = taxa[taxon]
                    raise ValueError(
                        "Duplicated taxon in accession file after alphanumeric cleanup: "
                        f"{taxon_raw} on line {row_idx}; first seen as {first_raw} on line {first_line}."
                    )
                taxa[taxon] = (row_idx, taxon_raw)

            if has_code_column:
                if len(row) < 2:
                    raise ValueError(f"Accession file line {row_idx} has no code column.")
                code = re.sub(r"\s+", "", row[1])
                if code:
                    if code in codes:
                        first_line, first_raw = codes[code]
                        raise ValueError(
                            f"Duplicated code in accession file: {code} on line {row_idx}; "
                            f"first seen for {first_raw} on line {first_line}."
                        )
                    codes[code] = (row_idx, taxon_raw)

    return taxa, codes


def validate_against_main(rows, main_clean_file: Optional[Path], has_code_column: bool) -> None:
    if main_clean_file is None:
        return

    main_taxa, main_codes = read_main_keys(main_clean_file, has_code_column)
    errors = []
    for row in rows:
        if row.taxon in main_taxa:
            main_line, main_raw = main_taxa[row.taxon]
            errors.append(
                f"Local taxon on line {row.line_no} duplicates accession-file taxon after "
                f"alphanumeric cleanup: {row.taxon_raw} -> {row.taxon}; "
                f"first seen as {main_raw} on accession line {main_line}."
            )
        if has_code_column and row.code in main_codes:
            main_line, main_raw = main_codes[row.code]
            errors.append(
                f"Local code on line {row.line_no} duplicates accession-file code: {row.code}; "
                f"first seen for {main_raw} on accession line {main_line}."
            )

    if errors:
        raise ValueError("\n".join(errors))


def validate_annotations(rows, feature_types, group_by) -> None:
    total_feature_rows = 0
    total_records = 0
    feature_type_counts = Counter()

    for row in rows:
        records, stats = make_seqrecords(row, feature_types, group_by, emit_warnings=False)
        if not records:
            raise ValueError(f"No local sequence units generated for local taxon {row.taxon_raw}")
        total_feature_rows += stats.feature_count
        total_records += len(records)
        feature_type_counts.update(stats.feature_type_counts)

    feature_summary = ", ".join(f"{feature}:{count}" for feature, count in sorted(feature_type_counts.items()))
    log_info(
        "Validated local assembly annotations: "
        f"taxa={len(rows)}, sequence_units={total_records}, feature_rows={total_feature_rows}, "
        f"feature_types={feature_summary or 'none'}, group_by={group_by or 'none'}"
    )


def write_local_taxa(rows, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as fh:
        for row in rows:
            fh.write(f"{row.taxon}\n")


def write_status(output: Path, has_code_column: bool, row_count: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"HAS_CODE_COLUMN={'true' if has_code_column else 'false'}\n")
        fh.write(f"LOCAL_ASSEMBLY_COUNT={row_count}\n")


def main() -> None:
    args = parse_args()
    try:
        feature_types = parse_feature_types(args.feature_types)
        group_by = parse_group_by(args.group_by)
        lines = read_manifest_lines(args.manifest)
        has_code_column = detect_has_code_column(lines, args.code_mode)
        rows = load_manifest(args.manifest, has_code_column)
        validate_against_main(rows, args.main_clean_file, has_code_column)
        validate_annotations(rows, feature_types, group_by)
        write_clean_manifest(rows, has_code_column, args.local_clean_output)
        write_local_taxa(rows, args.local_taxa_output)
        write_status(args.status_output, has_code_column, len(rows))
        log_info(f"Validated local assemblies manifest: {args.manifest}")
    except Exception as exc:
        log_error(str(exc))
        sys.exit(1)


if __name__ == "__main__":
    main()
