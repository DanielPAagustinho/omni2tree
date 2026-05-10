#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from Bio import SeqIO

from extract_local_cds_from_gff import (
    load_manifest as load_local_manifest,
    make_seqrecords,
    parse_feature_types,
    parse_group_by,
)


SCHEMA_VERSION = 1
DB_SUFFIX = "_cds_from_genomic.fna"


def clean_alnum(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", text or "")


def parse_bool(raw: str) -> bool:
    value = raw.lower()
    if value == "true":
        return True
    if value == "false":
        return False
    raise argparse.ArgumentTypeError(f"Expected true or false, got: {raw}")


def stable_digest(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def records_digest(records: Sequence[Any]) -> str:
    h = hashlib.sha256()
    for record in records:
        h.update(str(record.id).encode("utf-8"))
        h.update(b"\0")
        h.update(str(record.description).encode("utf-8"))
        h.update(b"\0")
        h.update(str(record.seq).encode("utf-8"))
        h.update(b"\0")
    return h.hexdigest()


def db_file_summary(path: Path) -> Optional[Dict[str, Any]]:
    if not path.is_file() or path.stat().st_size == 0:
        return None
    records = list(SeqIO.parse(str(path), "fasta"))
    return {
        "record_count": len(records),
        "records_sha256": records_digest(records),
    }


def read_accession_rows(path: Optional[Path], has_code_column: bool) -> List[Dict[str, Any]]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return []

    rows: List[Dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(line for line in fh if line.strip() and not line.lstrip().startswith("#"))
        try:
            next(reader)
        except StopIteration:
            return rows

        for row_idx, row in enumerate(reader, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            taxon_raw = row[0].strip()
            taxon = clean_alnum(taxon_raw)
            if not taxon:
                raise ValueError(f"Accession clean file line {row_idx} has an empty taxon after cleanup")
            if has_code_column:
                code = re.sub(r"\s+", "", row[1]) if len(row) > 1 else ""
                accession_cells = row[2:]
            else:
                code = None
                accession_cells = row[1:]
            accessions = [re.sub(r"\s+", "", value) for value in accession_cells if value.strip()]
            rows.append(
                {
                    "source": "ncbi",
                    "taxon": taxon,
                    "taxon_raw": taxon_raw,
                    "code": code,
                    "accessions": accessions,
                }
            )
    return rows


def read_local_rows(
    path: Optional[Path],
    has_code_column: bool,
    feature_types: Sequence[str],
    group_by: Optional[str],
) -> List[Dict[str, Any]]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return []

    rows = []
    for row in load_local_manifest(path, has_code_column):
        records, stats = make_seqrecords(row, feature_types, group_by)
        rows.append(
            {
                "source": "local",
                "taxon": row.taxon,
                "taxon_raw": row.taxon_raw,
                "code": row.code,
                "record_count": len(records),
                "records_sha256": records_digest(records),
                "feature_rows": stats.feature_count,
                "feature_type_counts": dict(sorted(stats.feature_type_counts.items())),
                "dna_path": str(row.dna_path),
                "annotation_path": str(row.annotation_path),
            }
        )
    return rows


def expected_taxa(ncbi_rows: Sequence[Dict[str, Any]], local_rows: Sequence[Dict[str, Any]]) -> Dict[str, str]:
    expected: Dict[str, str] = {}
    for row in [*ncbi_rows, *local_rows]:
        taxon = row["taxon"]
        if taxon in expected:
            raise ValueError(f"Taxon '{taxon}' appears in both current reference inputs")
        expected[taxon] = row["source"]
    return expected


def validate_no_extra_db_files(db_dir: Path, expected: Dict[str, str]) -> None:
    if not db_dir.is_dir():
        return
    extra = []
    for path in sorted(db_dir.glob(f"*{DB_SUFFIX}")):
        taxon = path.name[: -len(DB_SUFFIX)]
        if taxon not in expected:
            extra.append(path.name)
    if extra:
        raise ValueError(
            "Found file(s) in the 'db/' folder that are not present in the current reference inputs: "
            + ", ".join(extra)
            + ". Please move them away before using --resume."
        )


def db_outputs(db_dir: Path, expected: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
    outputs: Dict[str, Dict[str, Any]] = {}
    for taxon, source in sorted(expected.items()):
        filename = f"{taxon}{DB_SUFFIX}"
        path = db_dir / filename
        summary = db_file_summary(path)
        if summary is None:
            continue
        outputs[filename] = {
            "taxon": taxon,
            "source": source,
            **summary,
        }
    return outputs


def build_model(args: argparse.Namespace) -> Dict[str, Any]:
    feature_types = parse_feature_types(args.local_features)
    group_by = parse_group_by(args.local_group_by)
    ncbi_rows = read_accession_rows(args.clean_file, args.has_code_column) if args.input_enabled else []
    local_rows = (
        read_local_rows(args.local_clean_file, args.has_code_column, feature_types, group_by)
        if args.local_enabled
        else []
    )
    expected = expected_taxa(ncbi_rows, local_rows)
    if args.validate_db:
        validate_no_extra_db_files(args.db_dir, expected)

    ncbi_signature_payload = {
        "enabled": args.input_enabled,
        "params": {
            "use_mat_peptides": args.mat_peptides,
            "use_mat_peptides_only": args.only_mat_peptides,
        },
        "rows": sorted(
            (
                {
                    "taxon": row["taxon"],
                    "accessions": row["accessions"],
                }
                for row in ncbi_rows
            ),
            key=lambda item: item["taxon"],
        ),
    }
    local_signature_payload = {
        "enabled": args.local_enabled,
        "params": {
            "local_features": list(feature_types),
            "local_group_by": group_by,
        },
        "rows": sorted(
            (
                {
                    "taxon": row["taxon"],
                    "record_count": row["record_count"],
                    "records_sha256": row["records_sha256"],
                    "feature_rows": row["feature_rows"],
                    "feature_type_counts": row["feature_type_counts"],
                }
                for row in local_rows
            ),
            key=lambda item: item["taxon"],
        ),
    }
    if args.has_code_column:
        code_rows = sorted(
            (
                {"taxon": row["taxon"], "code": row["code"]}
                for row in [*ncbi_rows, *local_rows]
                if row.get("code")
            ),
            key=lambda item: item["taxon"],
        )
    else:
        code_rows = []

    codes_signature_payload = {
        "has_code_column": args.has_code_column,
        "rows": code_rows,
    }

    return {
        "schema_version": SCHEMA_VERSION,
        "complete": args.complete,
        "signatures": {
            "ncbi": stable_digest(ncbi_signature_payload),
            "local": stable_digest(local_signature_payload),
            "codes": stable_digest(codes_signature_payload),
        },
        "ncbi": {
            "enabled": args.input_enabled,
            "params": ncbi_signature_payload["params"],
            "rows": ncbi_rows,
        },
        "local": {
            "enabled": args.local_enabled,
            "params": local_signature_payload["params"],
            "rows": local_rows,
        },
        "codes": codes_signature_payload,
        "outputs": {
            "db": db_outputs(args.db_dir, expected),
        },
    }


def load_state(path: Path) -> Dict[str, Any]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Resume state file '{path}' is missing. Cannot use --resume.")
    with path.open(encoding="utf-8") as fh:
        state = json.load(fh)
    if state.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"Resume state file '{path}' has unsupported schema version: {state.get('schema_version')}"
        )
    if state.get("complete") is not True:
        raise ValueError(f"Resume state file '{path}' is incomplete. Cannot use --resume safely.")
    return state


def output_mismatch(
    previous: Dict[str, Any],
    current: Dict[str, Any],
    source: str,
) -> bool:
    previous_outputs = previous.get("outputs", {}).get("db", {})
    current_outputs = current.get("outputs", {}).get("db", {})
    current_expected_rows = current.get(source, {}).get("rows", [])
    for row in current_expected_rows:
        filename = f"{row['taxon']}{DB_SUFFIX}"
        old = previous_outputs.get(filename)
        new = current_outputs.get(filename)
        if old is None or new is None:
            return True
        if old.get("record_count") != new.get("record_count"):
            return True
        if old.get("records_sha256") != new.get("records_sha256"):
            return True
    return False


def bool_env(value: bool) -> str:
    return "true" if value else "false"


def cmd_check(args: argparse.Namespace) -> None:
    load_state(args.state)


def cmd_compare(args: argparse.Namespace) -> None:
    previous = load_state(args.state)
    current = build_model(args)

    ncbi_signature_changed = (
        previous.get("signatures", {}).get("ncbi") != current["signatures"]["ncbi"]
    )
    local_signature_changed = (
        previous.get("signatures", {}).get("local") != current["signatures"]["local"]
    )
    codes_changed = previous.get("signatures", {}).get("codes") != current["signatures"]["codes"]
    ncbi_db_mismatch = output_mismatch(previous, current, "ncbi")
    local_db_mismatch = output_mismatch(previous, current, "local")
    ncbi_rebuild = current["ncbi"]["enabled"] and (ncbi_signature_changed or ncbi_db_mismatch)
    local_rebuild = current["local"]["enabled"] and (local_signature_changed or local_db_mismatch)
    force_step4 = ncbi_rebuild or local_rebuild or codes_changed

    env = {
        "NCBI_SIGNATURE_CHANGED": ncbi_signature_changed,
        "LOCAL_SIGNATURE_CHANGED": local_signature_changed,
        "CODES_CHANGED": codes_changed,
        "NCBI_DB_MISMATCH": ncbi_db_mismatch,
        "LOCAL_DB_MISMATCH": local_db_mismatch,
        "NCBI_REBUILD": ncbi_rebuild,
        "LOCAL_REBUILD": local_rebuild,
        "FORCE_STEP4": force_step4,
    }
    with args.env_output.open("w", encoding="utf-8") as fh:
        for key, value in env.items():
            fh.write(f"{key}={bool_env(value)}\n")


def cmd_write(args: argparse.Namespace) -> None:
    state = build_model(args)
    args.state.parent.mkdir(parents=True, exist_ok=True)
    with args.state.open("w", encoding="utf-8") as fh:
        json.dump(state, fh, indent=2, sort_keys=True)
        fh.write("\n")


def cmd_write_codes(args: argparse.Namespace) -> None:
    if not args.has_code_column:
        return
    feature_types = parse_feature_types(args.local_features)
    group_by = parse_group_by(args.local_group_by)
    rows = []
    if args.input_enabled:
        rows.extend(read_accession_rows(args.clean_file, True))
    if args.local_enabled:
        rows.extend(read_local_rows(args.local_clean_file, True, feature_types, group_by))

    lines = []
    seen_taxa = set()
    seen_codes = set()
    for row in rows:
        taxon = row["taxon"]
        code = row.get("code")
        if not code:
            continue
        if taxon in seen_taxa:
            raise ValueError(f"Duplicated taxon while writing code mapping: {taxon}")
        if code in seen_codes:
            raise ValueError(f"Duplicated 5-letter code while writing code mapping: {code}")
        target = args.db_dir / f"{taxon}{DB_SUFFIX}"
        summary = db_file_summary(target)
        if summary is None or summary["record_count"] == 0:
            continue
        seen_taxa.add(taxon)
        seen_codes.add(code)
        lines.append(f"{taxon}\t{code}\n")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as fh:
        fh.writelines(lines)
    print(f"Wrote {len(lines)} current 5-letter code mapping(s) to {args.output}")


def add_model_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--db-dir", required=True, type=Path)
    parser.add_argument("--clean-file", type=Path, default=None)
    parser.add_argument("--local-clean-file", type=Path, default=None)
    parser.add_argument("--input-enabled", type=parse_bool, required=True)
    parser.add_argument("--local-enabled", type=parse_bool, required=True)
    parser.add_argument("--has-code-column", type=parse_bool, required=True)
    parser.add_argument("--mat-peptides", type=parse_bool, required=True)
    parser.add_argument("--only-mat-peptides", type=parse_bool, required=True)
    parser.add_argument("--local-features", default="CDS")
    parser.add_argument("--local-group-by", default=None)
    parser.add_argument("--validate-db", action="store_true")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Manage Omni2Tree step1 resume state.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    check = sub.add_parser("check", help="Validate that a resume state file exists and is complete")
    check.add_argument("--state", required=True, type=Path)
    check.set_defaults(func=cmd_check)

    compare = sub.add_parser("compare", help="Compare current inputs with a previous resume state")
    compare.add_argument("--state", required=True, type=Path)
    compare.add_argument("--env-output", required=True, type=Path)
    add_model_args(compare)
    compare.set_defaults(func=cmd_compare, complete=True)

    write = sub.add_parser("write", help="Write the current resume state")
    write.add_argument("--state", required=True, type=Path)
    write.add_argument("--complete", type=parse_bool, required=True)
    add_model_args(write)
    write.set_defaults(func=cmd_write)

    write_codes = sub.add_parser("write-codes", help="Write five_letter_taxon.tsv from current inputs and db files")
    write_codes.add_argument("--output", required=True, type=Path)
    add_model_args(write_codes)
    write_codes.set_defaults(func=cmd_write_codes, complete=True, validate_db=False)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        args.func(args)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
