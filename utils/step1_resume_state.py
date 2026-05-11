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
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
    except Exception:
        return None
    if not records:
        return None
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
        records, stats = make_seqrecords(row, feature_types, group_by, emit_warnings=False)
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


def stale_db_files(db_dir: Path, expected: Dict[str, str]) -> List[str]:
    if not db_dir.is_dir():
        return []
    extra = []
    for path in sorted(db_dir.glob(f"*{DB_SUFFIX}")):
        taxon = path.name[: -len(DB_SUFFIX)]
        if taxon not in expected:
            extra.append(path.name)
    return extra


def validate_no_extra_db_files(db_dir: Path, expected: Dict[str, str]) -> None:
    extra = stale_db_files(db_dir, expected)
    if extra:
        raise ValueError(
            "Found file(s) in the 'db/' folder that are not present in the current reference inputs: "
            + ", ".join(extra)
            + ". Please move them away before writing a complete resume state."
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


def validate_expected_db_outputs(outputs: Dict[str, Dict[str, Any]], expected: Dict[str, str]) -> None:
    missing = []
    for taxon in sorted(expected):
        filename = f"{taxon}{DB_SUFFIX}"
        if filename not in outputs:
            missing.append(filename)
    if missing:
        raise ValueError(
            "Missing, empty, or invalid expected db FASTA file(s): " + ", ".join(missing)
        )


def row_signatures(model: Dict[str, Any], source: str) -> Dict[str, str]:
    rows = model.get(source, {}).get("rows", [])
    signatures: Dict[str, str] = {}
    for row in rows:
        taxon = row["taxon"]
        if source == "ncbi":
            payload = {
                "taxon": taxon,
                "accessions": row.get("accessions", []),
            }
        elif source == "local":
            payload = {
                "taxon": taxon,
                "record_count": row.get("record_count"),
                "records_sha256": row.get("records_sha256"),
                "feature_rows": row.get("feature_rows"),
                "feature_type_counts": row.get("feature_type_counts", {}),
            }
        else:
            raise ValueError(f"Unsupported source for row signatures: {source}")
        signatures[taxon] = stable_digest(payload)
    return signatures


def source_params_signature(model: Dict[str, Any], source: str) -> str:
    payload = {
        "enabled": model.get(source, {}).get("enabled", False),
        "params": model.get(source, {}).get("params", {}),
    }
    return stable_digest(payload)


def output_by_taxon(model: Dict[str, Any], source: str) -> Dict[str, Dict[str, Any]]:
    outputs = model.get("outputs", {}).get("db", {})
    selected: Dict[str, Dict[str, Any]] = {}
    for item in outputs.values():
        if item.get("source") == source and item.get("taxon"):
            selected[item["taxon"]] = item
    return selected


def trusted_preserved_outputs(previous: Dict[str, Any], current: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    expected = expected_taxa(current["ncbi"]["rows"], current["local"]["rows"])
    preserved: Dict[str, Dict[str, Any]] = {}
    previous_rows = {
        "ncbi": row_signatures(previous, "ncbi"),
        "local": row_signatures(previous, "local"),
    }
    current_rows = {
        "ncbi": row_signatures(current, "ncbi"),
        "local": row_signatures(current, "local"),
    }
    params_match = {
        source: source_params_signature(previous, source) == source_params_signature(current, source)
        for source in ("ncbi", "local")
    }

    for filename, item in previous.get("outputs", {}).get("db", {}).items():
        taxon = item.get("taxon")
        source = item.get("source")
        if source not in {"ncbi", "local"} or not taxon:
            continue
        if expected.get(taxon) != source:
            continue
        if not params_match[source]:
            continue
        if previous_rows[source].get(taxon) != current_rows[source].get(taxon):
            continue
        if item.get("record_count") is None or item.get("records_sha256") is None:
            continue
        preserved[filename] = item
    return preserved


def repair_taxa(
    previous: Dict[str, Any],
    current: Dict[str, Any],
    source: str,
    params_changed: bool,
) -> Tuple[List[str], bool]:
    current_enabled = current.get(source, {}).get("enabled", False)
    if not current_enabled or params_changed:
        return [], False

    previous_rows = row_signatures(previous, source)
    current_rows = row_signatures(current, source)
    previous_outputs = output_by_taxon(previous, source)
    current_outputs = output_by_taxon(current, source)
    repair: List[str] = []

    for taxon in sorted(current_rows):
        current_output = current_outputs.get(taxon)
        old_signature = previous_rows.get(taxon)
        if old_signature is None:
            repair.append(taxon)
            continue
        if old_signature != current_rows[taxon]:
            repair.append(taxon)
            continue
        old_output = previous_outputs.get(taxon)
        if old_output is None:
            repair.append(taxon)
            continue
        if current_output is None:
            repair.append(taxon)
            continue
        if old_output.get("record_count") != current_output.get("record_count"):
            repair.append(taxon)
            continue
        if old_output.get("records_sha256") != current_output.get("records_sha256"):
            repair.append(taxon)
            continue

    rows_changed = previous_rows != current_rows
    return repair, rows_changed


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

    outputs = db_outputs(args.db_dir, expected)

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
            "db": outputs,
        },
    }


def load_state(path: Path, require_complete: bool = True) -> Dict[str, Any]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Resume state file '{path}' is missing. Cannot use --resume.")
    with path.open(encoding="utf-8") as fh:
        state = json.load(fh)
    if state.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"Resume state file '{path}' has unsupported schema version: {state.get('schema_version')}"
        )
    if require_complete and state.get("complete") is not True:
        raise ValueError(f"Resume state file '{path}' is incomplete. Cannot use --resume safely.")
    return state


def bool_env(value: bool) -> str:
    return "true" if value else "false"


def env_value(value: Any) -> str:
    if isinstance(value, bool):
        return bool_env(value)
    return str(value)


def cmd_check(args: argparse.Namespace) -> None:
    load_state(args.state, require_complete=False)


def cmd_compare(args: argparse.Namespace) -> None:
    previous = load_state(args.state, require_complete=False)
    current = build_model(args)

    ncbi_params_changed = (
        current["ncbi"]["enabled"]
        and previous.get("ncbi", {}).get("enabled")
        and source_params_signature(previous, "ncbi") != source_params_signature(current, "ncbi")
    )
    local_params_changed = (
        current["local"]["enabled"]
        and previous.get("local", {}).get("enabled")
        and source_params_signature(previous, "local") != source_params_signature(current, "local")
    )
    codes_changed = previous.get("signatures", {}).get("codes") != current["signatures"]["codes"]
    ncbi_rebuild = bool(ncbi_params_changed)
    local_rebuild = bool(local_params_changed)
    ncbi_repair, ncbi_rows_changed = repair_taxa(previous, current, "ncbi", ncbi_rebuild)
    local_repair, local_rows_changed = repair_taxa(previous, current, "local", local_rebuild)
    expected = expected_taxa(current["ncbi"]["rows"], current["local"]["rows"])
    stale_db = stale_db_files(args.db_dir, expected)
    ncbi_db_mismatch = bool(ncbi_repair)
    local_db_mismatch = bool(local_repair)
    force_step4 = (
        ncbi_rebuild
        or local_rebuild
        or codes_changed
        or bool(ncbi_repair)
        or bool(local_repair)
        or ncbi_rows_changed
        or local_rows_changed
        or bool(stale_db)
    )

    if args.ncbi_repair_output:
        args.ncbi_repair_output.write_text("".join(f"{taxon}\n" for taxon in ncbi_repair), encoding="utf-8")
    if args.local_repair_output:
        args.local_repair_output.write_text("".join(f"{taxon}\n" for taxon in local_repair), encoding="utf-8")
    if args.stale_db_output:
        args.stale_db_output.write_text("".join(f"{filename}\n" for filename in stale_db), encoding="utf-8")
    if args.expected_taxa_output:
        args.expected_taxa_output.write_text(
            "".join(f"{taxon}\t{source}\n" for taxon, source in sorted(expected.items())),
            encoding="utf-8",
        )

    env = {
        "NCBI_SIGNATURE_CHANGED": ncbi_params_changed or ncbi_rows_changed,
        "LOCAL_SIGNATURE_CHANGED": local_params_changed or local_rows_changed,
        "NCBI_PARAMS_CHANGED": ncbi_params_changed,
        "LOCAL_PARAMS_CHANGED": local_params_changed,
        "CODES_CHANGED": codes_changed,
        "NCBI_DB_MISMATCH": ncbi_db_mismatch,
        "LOCAL_DB_MISMATCH": local_db_mismatch,
        "NCBI_REBUILD": ncbi_rebuild,
        "LOCAL_REBUILD": local_rebuild,
        "NCBI_REPAIR_COUNT": len(ncbi_repair),
        "LOCAL_REPAIR_COUNT": len(local_repair),
        "NCBI_ROWS_CHANGED": ncbi_rows_changed,
        "LOCAL_ROWS_CHANGED": local_rows_changed,
        "STALE_DB_COUNT": len(stale_db),
        "FORCE_STEP4": force_step4,
    }
    with args.env_output.open("w", encoding="utf-8") as fh:
        for key, value in env.items():
            fh.write(f"{key}={env_value(value)}\n")


def cmd_write(args: argparse.Namespace) -> None:
    state = build_model(args)
    if args.complete:
        validate_expected_db_outputs(state["outputs"]["db"], expected_taxa(state["ncbi"]["rows"], state["local"]["rows"]))
    elif args.preserve_existing_outputs and args.state.is_file() and args.state.stat().st_size > 0:
        previous = load_state(args.state, require_complete=False)
        state["outputs"]["db"] = trusted_preserved_outputs(previous, state)
    args.state.parent.mkdir(parents=True, exist_ok=True)
    with args.state.open("w", encoding="utf-8") as fh:
        json.dump(state, fh, indent=2, sort_keys=True)
        fh.write("\n")


def cmd_mark_output(args: argparse.Namespace) -> None:
    current = build_model(args)
    expected = expected_taxa(current["ncbi"]["rows"], current["local"]["rows"])
    taxon = clean_alnum(args.taxon)
    if not taxon:
        raise ValueError("Cannot mark output for an empty taxon")
    if expected.get(taxon) != args.source:
        raise ValueError(
            f"Taxon '{taxon}' is not a current {args.source} reference input; cannot mark db output."
        )

    target = args.db_dir / f"{taxon}{DB_SUFFIX}"
    summary = db_file_summary(target)
    if summary is None:
        raise ValueError(f"Cannot mark missing, empty, or invalid db FASTA output: {target}")

    if args.state.is_file() and args.state.stat().st_size > 0:
        previous = load_state(args.state, require_complete=False)
        outputs = trusted_preserved_outputs(previous, current)
    else:
        outputs = {}

    filename = f"{taxon}{DB_SUFFIX}"
    outputs[filename] = {
        "taxon": taxon,
        "source": args.source,
        **summary,
    }
    current["complete"] = False
    current["outputs"]["db"] = outputs
    args.state.parent.mkdir(parents=True, exist_ok=True)
    with args.state.open("w", encoding="utf-8") as fh:
        json.dump(current, fh, indent=2, sort_keys=True)
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
    compare.add_argument("--ncbi-repair-output", type=Path, default=None)
    compare.add_argument("--local-repair-output", type=Path, default=None)
    compare.add_argument("--stale-db-output", type=Path, default=None)
    compare.add_argument("--expected-taxa-output", type=Path, default=None)
    add_model_args(compare)
    compare.set_defaults(func=cmd_compare, complete=True)

    write = sub.add_parser("write", help="Write the current resume state")
    write.add_argument("--state", required=True, type=Path)
    write.add_argument("--complete", type=parse_bool, required=True)
    write.add_argument("--preserve-existing-outputs", action="store_true")
    add_model_args(write)
    write.set_defaults(func=cmd_write)

    mark_output = sub.add_parser("mark-output", help="Mark one completed db FASTA output in the resume state")
    mark_output.add_argument("--state", required=True, type=Path)
    mark_output.add_argument("--taxon", required=True)
    mark_output.add_argument("--source", choices=("ncbi", "local"), required=True)
    add_model_args(mark_output)
    mark_output.set_defaults(func=cmd_mark_output, complete=False, validate_db=False)

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
