#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import unquote

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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


def clean_alnum(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", text or "")


def clean_token(text: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]", "_", text or "")
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "unknown"


MISSING_VALUES = {"", ".", "-", "na", "n/a", "unknown"}


def is_missing_value(value: Optional[str]) -> bool:
    return value is None or value.strip().strip('"').strip("'").lower() in MISSING_VALUES


def normalize_attr_value(value: Optional[str]) -> str:
    if is_missing_value(value):
        return ""
    return value.strip().strip('"').strip("'")


def parse_feature_types(raw: str) -> Tuple[str, ...]:
    values = tuple(part.strip() for part in raw.split(",") if part.strip())
    if not values:
        raise ValueError("At least one local feature type must be provided.")
    if any(re.search(r"\s", value) for value in values):
        raise ValueError("Local feature types must be comma-separated values without spaces inside each value.")
    return values


def parse_group_by(raw: Optional[str]) -> Optional[str]:
    if raw is None or is_missing_value(raw):
        return None
    value = raw.strip()
    if "," in value or re.search(r"\s", value):
        raise ValueError("Only one --group-by attribute is allowed.")
    return value


@dataclass
class ManifestRow:
    line_no: int
    taxon_raw: str
    taxon: str
    code: Optional[str]
    dna_path: Path
    annotation_path: Path


@dataclass
class CdsFeature:
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    phase: Optional[int]
    attrs: Dict[str, str]
    line_no: int


@dataclass
class ExtractionStats:
    feature_count: int
    record_count: int
    feature_type_counts: Dict[str, int]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract selected local GTF/GFF3 features into Omni2Tree db FASTA files."
    )
    p.add_argument("--manifest", required=True, type=Path, help="CSV with local assemblies")
    p.add_argument("--db-dir", required=True, type=Path, help="Output db directory")
    p.add_argument(
        "--has-code-column",
        action="store_true",
        help="Manifest is expected to contain taxon,code,dna,gff columns",
    )
    p.add_argument(
        "--codes-output",
        type=Path,
        default=None,
        help="five_letter_taxon.tsv to update when --has-code-column is used",
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="Skip local taxa whose db FASTA already exists and is non-empty",
    )
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
    return p.parse_args()


def resolve_path(raw_path: str, base_dir: Path) -> Path:
    path = Path(raw_path.strip())
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def load_manifest(path: Path, has_code_column: bool) -> List[ManifestRow]:
    base_dir = path.resolve().parent
    rows: List[ManifestRow] = []

    with path.open(newline="", encoding="utf-8") as fh:
        filtered = (line for line in fh if line.strip() and not line.lstrip().startswith("#"))
        reader = csv.reader(filtered)
        try:
            next(reader)
        except StopIteration as exc:
            raise ValueError(f"Local assemblies manifest is empty: {path}") from exc

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

            if has_code_column:
                code = row[1].strip()
                dna_raw = row[2]
                ann_raw = row[3]
            else:
                code = None
                dna_raw = row[1]
                ann_raw = row[2]

            dna_path = resolve_path(dna_raw, base_dir)
            annotation_path = resolve_path(ann_raw, base_dir)

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


def parse_gff3_attributes(field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if is_missing_value(field):
        return attrs
    for item in field.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            key, value = item, ""
        key = key.strip()
        if is_missing_value(key):
            continue
        values = [normalize_attr_value(unquote(v.strip())) for v in value.split(",")]
        values = [v for v in values if v]
        if values:
            attrs[key] = values[0]
    return attrs


def parse_gtf_attributes(field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if is_missing_value(field):
        return attrs
    for item in field.split(";"):
        item = item.strip()
        if not item:
            continue
        match = re.match(r"(\S+)\s+\"?([^\"]*)\"?$", item)
        if match:
            key = match.group(1).strip()
            value = normalize_attr_value(match.group(2))
            if not is_missing_value(key) and value:
                attrs[key] = value
        elif "=" in item:
            key, value = item.split("=", 1)
            key = key.strip()
            value = normalize_attr_value(value)
            if not is_missing_value(key) and value:
                attrs[key] = value
    return attrs


def parse_attributes(field: str) -> Dict[str, str]:
    gff3_attrs = parse_gff3_attributes(field)
    gtf_attrs = parse_gtf_attributes(field)
    attrs = dict(gff3_attrs)
    attrs.update({k: v for k, v in gtf_attrs.items() if v})
    return attrs


def read_features(path: Path, feature_types: Sequence[str]) -> List[CdsFeature]:
    wanted = {feature_type.lower() for feature_type in feature_types}
    features: List[CdsFeature] = []
    with path.open(encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                log_warn(f"Skipping malformed annotation line {line_no} in {path}: expected at least 9 columns")
                continue
            seqid, _source, feature_type, start, end, _score, strand, phase = cols[:8]
            attributes = "\t".join(cols[8:])
            if feature_type.lower() not in wanted:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                log_warn(f"Skipping feature line {line_no} in {path}: invalid coordinates")
                continue
            if start_i < 1 or end_i < start_i:
                log_warn(f"Skipping feature line {line_no} in {path}: unsupported coordinate range {start}..{end}")
                continue
            if strand not in {"+", "-"}:
                log_warn(f"Feature line {line_no} in {path} has strand '{strand}', treating as '+'")
                strand = "+"
            phase_i: Optional[int]
            if phase in {".", ""}:
                phase_i = None
            else:
                try:
                    phase_i = int(phase)
                except ValueError:
                    phase_i = None
                if phase_i not in {0, 1, 2}:
                    phase_i = None

            features.append(
                CdsFeature(
                    seqid=seqid,
                    feature_type=feature_type,
                    start=start_i,
                    end=end_i,
                    strand=strand,
                    phase=phase_i,
                    attrs=parse_attributes(attributes),
                    line_no=line_no,
                )
            )
    return features


def attr_value(attrs: Dict[str, str], keys: Iterable[str]) -> str:
    for key in keys:
        value = normalize_attr_value(attrs.get(key))
        if value:
            return value
    return ""


def group_key(feature: CdsFeature, group_by: Optional[str]) -> str:
    if group_by:
        value = attr_value(feature.attrs, (group_by,))
        if not value:
            raise ValueError(
                f"Feature line {feature.line_no} is missing group-by attribute '{group_by}' "
                "or its value is missing."
            )
        return value
    return f"{feature.seqid}:{feature.start}-{feature.end}:{feature.strand}:line{feature.line_no}"


def load_fasta(path: Path) -> Dict[str, SeqRecord]:
    records = SeqIO.to_dict(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def build_location(features: List[CdsFeature]) -> str:
    parts = [f"{feat.seqid}:{feat.start}..{feat.end}({feat.strand})" for feat in features]
    return ",".join(parts)


def extract_group_sequence(features: List[CdsFeature], fasta_records: Dict[str, SeqRecord]):
    strand = features[0].strand
    ordered = sorted(features, key=lambda feat: feat.start, reverse=(strand == "-"))
    parts = []
    for feat in ordered:
        record = fasta_records.get(feat.seqid)
        if record is None:
            raise ValueError(f"Sequence ID '{feat.seqid}' from annotation line {feat.line_no} not found in FASTA")
        if feat.end > len(record.seq):
            raise ValueError(
                f"Feature coordinates {feat.seqid}:{feat.start}..{feat.end} exceed FASTA sequence length "
                f"({len(record.seq)})"
            )
        part = record.seq[feat.start - 1 : feat.end]
        if strand == "-":
            part = part.reverse_complement()
        if feat.phase in {1, 2}:
            part = part[feat.phase :]
        parts.append(part)
    return sum(parts[1:], parts[0]) if parts else ""


def make_seqrecords(
    row: ManifestRow,
    feature_types: Sequence[str] = ("CDS",),
    group_by: Optional[str] = None,
) -> Tuple[List[SeqRecord], ExtractionStats]:
    fasta_records = load_fasta(row.dna_path)
    features = read_features(row.annotation_path, feature_types)
    if not features:
        raise ValueError(
            f"No selected feature(s) found in {row.annotation_path} for local taxon {row.taxon_raw}: "
            f"{', '.join(feature_types)}"
        )

    grouped: Dict[str, List[CdsFeature]] = {}
    for feat in features:
        grouped.setdefault(group_key(feat, group_by), []).append(feat)

    seqrecords: List[SeqRecord] = []
    used_protein_ids: Dict[str, int] = {}
    source_name = clean_token(row.dna_path.stem)

    for idx, key in enumerate(sorted(grouped), start=1):
        group = grouped[key]
        first = group[0]
        seq = extract_group_sequence(group, fasta_records)
        if not seq:
            raise ValueError(f"Empty local feature sequence generated for local taxon {row.taxon_raw}, group {key}")
        if len(seq) % 3 != 0:
            log_warn(
                f"Local feature group '{key}' for taxon {row.taxon_raw} has length {len(seq)}, "
                "which is not divisible by 3. Step 1 cleaning will pad it with N."
            )

        gene = attr_value(first.attrs, ("gene", "gene_name", "Name", "locus_tag", "gene_id"))
        product = attr_value(first.attrs, ("product", "protein", "Name"))
        locus_tag = attr_value(first.attrs, ("locus_tag", "gene_id"))
        db_xref = attr_value(first.attrs, ("Dbxref", "db_xref"))
        protein_id = attr_value(first.attrs, ("protein_id", "proteinId"))
        if not protein_id and group_by:
            protein_id = attr_value(first.attrs, (group_by,))
        if not protein_id:
            protein_id = f"{row.taxon}_local_cds_{idx}"

        if protein_id in used_protein_ids:
            used_protein_ids[protein_id] += 1
            protein_id_unique = f"{protein_id}_{used_protein_ids[protein_id]}"
        else:
            used_protein_ids[protein_id] = 1
            protein_id_unique = protein_id

        header = f"lcl|{source_name}_cds_{clean_token(protein_id_unique)}_{idx}"
        description_parts = []
        if gene:
            description_parts.append(f"[gene={gene}]")
        if product:
            description_parts.append(f"[protein={product}]")
        description_parts.append(f"[protein_id={protein_id_unique}]")
        description_parts.append(f"[location={build_location(group)}]")
        if locus_tag:
            description_parts.append(f"[locus_tag={locus_tag}]")
        if db_xref:
            description_parts.append(f"[db_xref={db_xref}]")
        description_parts.append(f"[gbkey={first.feature_type}]")

        seqrecords.append(
            SeqRecord(
                seq=seq,
                id=header,
                name=header,
                description=" ".join(description_parts),
            )
        )

    stats = ExtractionStats(
        feature_count=len(features),
        record_count=len(seqrecords),
        feature_type_counts=dict(Counter(feat.feature_type for feat in features)),
    )
    return seqrecords, stats


def load_codes(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not path.exists():
        return mapping
    with path.open(encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Invalid existing code mapping line {line_no}: {line}")
            mapping[parts[0]] = parts[1]
    return mapping


def prepare_code_updates(rows: List[ManifestRow], codes_output: Optional[Path]) -> List[tuple[str, str]]:
    if codes_output is None:
        return []
    existing = load_codes(codes_output)
    existing_codes = {code: taxon for taxon, code in existing.items()}
    to_append = []

    for row in rows:
        if row.code is None:
            continue
        if row.taxon in existing:
            if existing[row.taxon] != row.code:
                raise ValueError(
                    f"Taxon {row.taxon_raw} already has code {existing[row.taxon]} in {codes_output}, "
                    f"but local manifest provides {row.code}"
                )
            continue
        if row.code in existing_codes:
            raise ValueError(
                f"Local code {row.code} for taxon {row.taxon_raw} is already used by "
                f"{existing_codes[row.code]} in {codes_output}"
            )
        to_append.append((row.taxon, row.code))

    return to_append


def append_codes(codes_output: Optional[Path], to_append: List[tuple[str, str]]) -> None:
    if codes_output is None or not to_append:
        return

    codes_output.parent.mkdir(parents=True, exist_ok=True)
    with codes_output.open("a", encoding="utf-8") as fh:
        for taxon, code in to_append:
            fh.write(f"{taxon}\t{code}\n")
    if to_append:
        log_info(f"Updated {codes_output} with {len(to_append)} local taxon code(s)")


def write_local_db_files(
    rows: List[ManifestRow],
    db_dir: Path,
    resume: bool,
    feature_types: Sequence[str],
    group_by: Optional[str],
) -> int:
    db_dir.mkdir(parents=True, exist_ok=True)
    written = 0
    skipped = 0
    total_feature_rows = 0
    total_records = 0
    feature_type_counts = Counter()
    for row in rows:
        target = db_dir / f"{row.taxon}_cds_from_genomic.fna"
        if resume and target.exists() and target.stat().st_size > 0:
            log_info(f"Local taxon {row.taxon_raw} already exists in db; skipping extraction: {target}")
            skipped += 1
            continue
        records, stats = make_seqrecords(row, feature_types, group_by)
        if not records:
            raise ValueError(f"No local sequence units generated for local taxon {row.taxon_raw}")
        with target.open("w", encoding="utf-8") as fh:
            SeqIO.write(records, fh, "fasta")
        written += 1
        total_feature_rows += stats.feature_count
        total_records += stats.record_count
        feature_type_counts.update(stats.feature_type_counts)
        log_info(
            f"Wrote {len(records)} local sequence unit(s) from {stats.feature_count} feature row(s) "
            f"for {row.taxon_raw}: {target}"
        )
    feature_summary = ", ".join(f"{feature}:{count}" for feature, count in sorted(feature_type_counts.items()))
    log_info(
        "Local extraction report: "
        f"taxa_written={written}, taxa_skipped={skipped}, feature_rows={total_feature_rows}, "
        f"sequence_units={total_records}, feature_types={feature_summary or 'none'}, "
        f"group_by={group_by or 'none'}"
    )
    log_info(f"Local assembly extraction summary: {written} written, {skipped} skipped")
    return written


def main() -> None:
    args = parse_args()
    try:
        feature_types = parse_feature_types(args.feature_types)
        group_by = parse_group_by(args.group_by)
        rows = load_manifest(args.manifest, args.has_code_column)
        log_info(f"Loaded {len(rows)} local assembly row(s) from {args.manifest}")
        code_updates = prepare_code_updates(rows, args.codes_output)
        written = write_local_db_files(rows, args.db_dir, args.resume, feature_types, group_by)
        append_codes(args.codes_output, code_updates)
        log_info(f"Local feature extraction completed successfully. Local db files written: {written}")
    except Exception as exc:
        log_error(str(exc))
        sys.exit(1)


if __name__ == "__main__":
    main()
