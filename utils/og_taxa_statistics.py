#!/usr/bin/env python
import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def og_sort_key(value):
    if value.startswith("OG") and value[2:].isdigit():
        return (0, int(value[2:]))
    return (1, value)


def main():
    p = argparse.ArgumentParser(
        description="Generate per-OG taxa counts table and histogram from an OG summary table."
    )
    p.add_argument("--input", required=True,
                   help="Input TSV/CSV with OG and Taxa columns")
    p.add_argument("--out-dir", default=".",
                   help="Output directory")
    p.add_argument("--prefix", default="og_taxa_per_og",
                   help="Prefix for the output files (without extension)")
    p.add_argument("--og-column", default="OG",
                   help="Column containing the orthologous group identifier")
    p.add_argument("--taxa-column", default="Taxa",
                   help="Column containing the taxa list")
    p.add_argument("--delimiter", default=";",
                   help="Delimiter used inside the taxa column")
    args = p.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    sep = "\t" if input_path.suffix.lower() in {".tsv", ".tab"} else ","
    df = pd.read_csv(input_path, sep=sep, dtype=str).fillna("")

    if args.og_column not in df.columns:
        raise SystemExit(f"Missing OG column '{args.og_column}' in {input_path}")
    if args.taxa_column not in df.columns:
        raise SystemExit(f"Missing taxa column '{args.taxa_column}' in {input_path}")

    taxa_by_og = defaultdict(set)
    for _, row in df.iterrows():
        og = str(row[args.og_column]).strip()
        taxa_field = str(row[args.taxa_column]).strip()
        if not og or not taxa_field:
            continue
        for taxon in taxa_field.split(args.delimiter):
            taxon = taxon.strip()
            if taxon:
                taxa_by_og[og].add(taxon)

    if not taxa_by_og:
        raise SystemExit(f"No OG/taxa data found in {input_path}")

    rows = []
    for og in sorted(taxa_by_og, key=og_sort_key):
        taxa = sorted(taxa_by_og[og])
        rows.append({
            "OG": og,
            "taxa_count": len(taxa),
            "taxa": args.delimiter.join(taxa),
        })

    counts_df = pd.DataFrame(rows)
    counts_df.to_csv(out_dir / f"{args.prefix}.tsv", sep="\t", index=False)

    data = counts_df["taxa_count"]
    UNIQUE_THRESHOLD = 10
    MAX_BINS = 30

    unique_counts = np.sort(data.unique())
    n_unique = len(unique_counts)

    fig, ax = plt.subplots(figsize=(10, 6))

    if n_unique <= UNIQUE_THRESHOLD:
        freq = data.value_counts().sort_index()
        ax.bar(freq.index, freq.values)
        ax.set_xticks(freq.index)
    else:
        bin_edges = np.histogram_bin_edges(data, bins="auto")
        if len(bin_edges) - 1 > MAX_BINS:
            bin_edges = np.linspace(data.min(), data.max(), MAX_BINS + 1)
        ax.hist(data, bins=bin_edges)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=15))

    ax.set_xlabel("Number of taxa")
    ax.set_ylabel("Number of orthologous groups")
    ax.set_title("Distribution of Taxa Counts per Orthologous Group")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.prefix}_distribution.png", dpi=300)

    freq_table = data.value_counts().sort_index()
    freq_df = freq_table.to_frame(name="n_ogs")
    freq_df["percent"] = (freq_df["n_ogs"] / freq_df["n_ogs"].sum() * 100).round(2)
    freq_df.index.name = "taxa_count"
    freq_df.to_csv(out_dir / f"{args.prefix}_frequency.tsv", sep="\t", header=True)


if __name__ == "__main__":
    main()
