#!/usr/bin/env python3
"""
Build an OG-to-gene mapping table for entropy analysis step 1 from Read2Tree OG_genes.tsv.

Rules for selecting one gene name per OG:
1. Choose the most frequent gene name within the OG.
2. If there is a tie, choose the alphabetically first gene name.

Output is a CSV with columns: OG,gene
This matches the format expected by msa_to_position_table.py in the entropy pipeline.
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare entropy step 1 OG->gene mapping from OG_genes.tsv"
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Read2Tree OG_genes table (TSV/CSV), e.g. OG_genes.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output CSV for entropy step 1 mapping (default: "
            "<input_stem>_entropy_step1_og_gene_table.csv)"
        ),
    )
    parser.add_argument(
        "--keep_empty_gene",
        action="store_true",
        help="Keep rows with empty gene names (default: drop empty gene values)",
    )
    return parser.parse_args()


def detect_separator(path: Path) -> str | None:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return "\t"
    if suffix == ".csv":
        return ","
    return None


def load_og_genes_table(path: Path) -> pd.DataFrame:
    sep = detect_separator(path)
    if sep is None:
        df = pd.read_csv(path, sep=None, engine="python", keep_default_na=False)
    else:
        df = pd.read_csv(path, sep=sep, keep_default_na=False)

    col_map = {col.lower(): col for col in df.columns}
    if "og" not in col_map:
        raise ValueError("Input table does not contain an 'OG' column")
    if "gene" not in col_map:
        raise ValueError("Input table does not contain a 'Gene' column")

    out = df[[col_map["og"], col_map["gene"]]].copy()
    out.columns = ["OG", "Gene"]
    out["OG"] = out["OG"].astype("string").str.strip()
    out["Gene"] = out["Gene"].astype("string").str.strip()
    return out


def og_sort_key(og_value: str):
    match = re.match(r"^([A-Za-z_]+)(\d+)$", og_value)
    if match:
        return (match.group(1), int(match.group(2)))
    return (og_value, float("inf"))


def choose_gene_name(gene_series: pd.Series) -> tuple[str, int, int]:
    counts = gene_series.value_counts(dropna=False)
    max_count = int(counts.max())
    winners = sorted([str(g) for g, c in counts.items() if int(c) == max_count], key=str.casefold)
    chosen = winners[0]
    return chosen, max_count, len(counts)


def build_mapping(df: pd.DataFrame, keep_empty_gene: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = df.copy()
    if not keep_empty_gene:
        work = work[work["Gene"].notna()]
        work = work[work["Gene"] != ""]

    if work.empty:
        raise ValueError("No valid OG/Gene rows available after filtering")

    records = []
    conflict_records = []

    for og, group in work.groupby("OG", sort=False):
        chosen, chosen_count, n_unique = choose_gene_name(group["Gene"])
        records.append({"OG": og, "gene": chosen})

        if n_unique > 1:
            gene_counts = group["Gene"].value_counts()
            for gene_name, count in gene_counts.items():
                conflict_records.append(
                    {
                        "OG": og,
                        "gene_name": str(gene_name),
                        "count": int(count),
                        "chosen_gene": chosen,
                        "chosen_count": chosen_count,
                    }
                )

    mapping_df = pd.DataFrame(records).sort_values("OG", key=lambda s: s.map(og_sort_key)).reset_index(drop=True)
    conflicts_df = pd.DataFrame(conflict_records)
    return mapping_df, conflicts_df


def main():
    args = parse_args()

    if not args.input.exists():
        log_error(f"Input file not found: {args.input}")
        sys.exit(1)

    output_path = args.output
    if output_path is None:
        output_path = args.input.with_name(f"{args.input.stem}_entropy_step1_og_gene_table.csv")

    try:
        og_df = load_og_genes_table(args.input)
        ambiguous_mask = og_df["Gene"].notna() & og_df["Gene"].str.upper().isin(["NA", "NAN"])
        if ambiguous_mask.any():
            ambiguous_ogs = sorted(og_df.loc[ambiguous_mask, "OG"].dropna().unique())
            log_warn(
                "Found gene names equal to NA/NAN; treating them as valid gene names "
                f"(OGs: {', '.join(ambiguous_ogs)})"
            )

        if not args.keep_empty_gene:
            empty_mask = og_df["Gene"].isna() | (og_df["Gene"] == "")
            if empty_mask.any():
                log_warn(f"Skipping {int(empty_mask.sum())} row(s) with empty gene name.")

        mapping_df, conflicts_df = build_mapping(og_df, args.keep_empty_gene)
    except Exception as exc:
        log_error(str(exc))
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_df.to_csv(output_path, index=False)

    log_info(f"Saved entropy step 1 OG->gene table: {output_path}")
    log_info(f"OG rows written: {len(mapping_df)}")
    log_info(f"Unique OGs in input: {og_df['OG'].nunique()}")

    if conflicts_df.empty:
        log_info("No OGs with multiple gene names detected.")
    else:
        n_conflict_ogs = conflicts_df["OG"].nunique()
        log_info(f"Detected {n_conflict_ogs} OG(s) with multiple gene names.")
        log_info("Selection rule applied: most frequent gene name; alphabetical tie-break.")


if __name__ == "__main__":
    main()
