#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def main():
    p = argparse.ArgumentParser(
        description="Generates TSV and plot of number of sequence units per reference taxon."
    )
    p.add_argument("--db-dir", required=True,
                   help="Directory that contains the *_cds_from_genomic.fna")
    p.add_argument("--out-dir", default=".",
                   help="Output directory")
    p.add_argument("--prefix", default="cds_count_per_accession",
                   help="Prefix for the output files (without extension)")
    args = p.parse_args()

    db_dir = Path(args.db_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(db_dir.glob("*_cds_from_genomic.fna"))
    if not files:
        raise SystemExit(f"No FASTA files found in {db_dir}")

    rows = []
    for f in files:
        # Counts FASTA headers
        with f.open() as fh:
            c = sum(1 for line in fh if line.startswith(">"))
        accession = f.name.replace("_cds_from_genomic.fna", "")
        rows.append((c, accession))

    df = pd.DataFrame(rows, columns=["cds_count", "accession"])
    tsv_path = out_dir / f"{args.prefix}.tsv"
    df.sort_values("cds_count").to_csv(tsv_path, sep="\t", header=False, index=False)

    data = df["cds_count"]
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
        bin_edges = np.histogram_bin_edges(data, bins='auto')
        if len(bin_edges) - 1 > MAX_BINS:
            bin_edges = np.linspace(data.min(), data.max(), MAX_BINS + 1)
        ax.hist(data, bins=bin_edges)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=15))

    ax.set_xlabel("Number of sequence units")
    ax.set_ylabel("Number of reference taxa")
    ax.set_title("Distribution of Sequence-Unit Counts per Reference Taxon")
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"{args.prefix}_distribution.png", dpi=300)

    # Frequency tables
    freq_table = data.value_counts().sort_index()
    freq_df = freq_table.to_frame(name="n_reference_taxa")
    freq_df["percent"] = (freq_df["n_reference_taxa"] / freq_df["n_reference_taxa"].sum() * 100).round(2)
    freq_df.index.name = "cds_count"
    freq_df.to_csv(out_dir / f"{args.prefix}_frequency.tsv", sep="\t", header=True)

if __name__ == "__main__":
    main()
