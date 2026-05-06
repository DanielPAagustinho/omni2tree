#!/usr/bin/env python3
"""
Calculate Shannon entropy from MSA position table.
Takes output from msa_to_position_table.py and calculates entropy per position.

Usage:
    python calculate_entropy.py --input positions.csv --output entropy.csv
    python calculate_entropy.py --input positions.csv --output entropy.csv --metadata meta.csv --group_by genotype

Author: Daniel Agustinho
Date: 2025
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys

def calculate_shannon_entropy(char_series):
    """
    Calculate Shannon entropy for a series of characters.
    
    H = -sum(p_i * log2(p_i)) where p_i is probability of character i
    """
    # Count frequencies
    counts = char_series.value_counts()
    total = len(char_series)
    
    if total == 0:
        return 0.0
    
    # Calculate probabilities
    probs = counts / total
    
    # Calculate entropy (in bits)
    entropy = -np.sum(probs * np.log2(probs))
    
    return entropy

def process_entropy(position_df):
    """Calculate entropy and statistics for a position group."""
    
    entropy = calculate_shannon_entropy(position_df['character'])
    n_samples = len(position_df)
    n_unique = position_df['character'].nunique()
    
    # Calculate gap percentage
    gap_count = (position_df['character'] == '-').sum()
    gap_pct = (gap_count / n_samples * 100) if n_samples > 0 else 0
    
    # Most common character
    most_common = position_df['character'].mode()[0] if len(position_df) > 0 else None
    
    return pd.Series({
        'entropy': entropy,
        'n_samples': n_samples,
        'n_unique_chars': n_unique,
        'gap_percent': gap_pct,
        'most_common_char': most_common
    })


def select_position_column(pos_df, use_reference_positions):
    """Choose the coordinate column for entropy grouping."""
    if use_reference_positions:
        if 'position_reference' not in pos_df.columns:
            print("\nError: --use_reference_positions specified but position_reference column not found", file=sys.stderr)
            print("       Run msa_to_position_table.py with --reference_id to generate reference positions.", file=sys.stderr)
            sys.exit(1)
        print("\n*** Using reference-based positions (position_reference column) ***")
        print("This provides ungapped positions in the selected reference sequence")
        return 'position_reference'

    if 'position' not in pos_df.columns:
        print("Error: Position table does not contain required 'position' column", file=sys.stderr)
        sys.exit(1)

    if 'position_reference' in pos_df.columns:
        print("\n*** Using alignment positions (position column) ***")
        print("  Reference positions are available; pass --use_reference_positions to use them")

    return 'position'

def main():
    parser = argparse.ArgumentParser(
        description='Calculate Shannon entropy from position table',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--input', required=True, type=Path,
                        help='Position table CSV from msa_to_position_table.py')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output entropy CSV file')
    parser.add_argument('--metadata', type=Path,
                        help='Metadata CSV for grouping')
    parser.add_argument('--group_by', nargs='+',
                        help='Columns to group by (must be in metadata)')
    parser.add_argument('--min_samples', type=int, default=5,
                        help='Minimum samples per position (default: 5)')
    parser.add_argument('--use_reference_positions', action='store_true',
                        help='Use position_reference instead of alignment position. Requires msa_to_position_table.py --reference_id.')
    parser.add_argument('--exclude_gaps', action='store_true',
                        help='Exclude gaps (-) from entropy calculation')
    
    args = parser.parse_args()
    
    # Check input exists
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Load position table
    print(f"Loading position table: {args.input}")
    pos_df = pd.read_csv(args.input)
    
    print(f"  Records: {len(pos_df):,}")
    print(f"  Samples: {pos_df['label'].nunique()}")
    print(f"  Genes: {pos_df['gene'].nunique()}")

    position_col = select_position_column(pos_df, args.use_reference_positions)

    if args.use_reference_positions:
        before = len(pos_df)
        pos_df = pos_df[pos_df[position_col].notna()]
        filtered = before - len(pos_df)
        if filtered > 0:
            print(f"  Filtered {filtered:,} records at gap positions in the reference sequence")
        print(f"  Remaining records: {len(pos_df):,}")
        if pos_df.empty:
            print("Error: No records remain after filtering reference-gap positions", file=sys.stderr)
            sys.exit(1)
    
    # Load and merge metadata if provided
    if args.metadata:
        if not args.metadata.exists():
            print(f"Error: Metadata not found: {args.metadata}", file=sys.stderr)
            sys.exit(1)
        
        print(f"Loading metadata: {args.metadata}")
        meta_df = pd.read_csv(args.metadata)
        pos_df = pos_df.merge(meta_df, on='label', how='left')
        
        # Check for grouping columns
        if args.group_by:
            for col in args.group_by:
                if col not in pos_df.columns:
                    print(f"Error: Column '{col}' not in data", file=sys.stderr)
                    sys.exit(1)
    
    # Exclude gaps if requested
    if args.exclude_gaps:
        before = len(pos_df)
        pos_df = pos_df[pos_df['character'] != '-']
        print(f"Excluded {before - len(pos_df):,} gaps")
    
    # Set up grouping columns
    group_cols = ['gene', 'og', position_col, 'seq_type']
    if args.group_by:
        group_cols = args.group_by + group_cols
    
    print(f"\nCalculating entropy...")
    print(f"  Grouping by: {', '.join(group_cols)}")
    
    # Calculate entropy
    entropy_df = pos_df.groupby(group_cols).apply(process_entropy).reset_index()
    if position_col != 'position':
        entropy_df = entropy_df.rename(columns={position_col: 'position'})
    
    # Filter by minimum samples
    before = len(entropy_df)
    entropy_df = entropy_df[entropy_df['n_samples'] >= args.min_samples]
    filtered = before - len(entropy_df)
    if filtered > 0:
        print(f"  Filtered {filtered} positions with < {args.min_samples} samples")
    if entropy_df.empty:
        print("Error: Entropy table is empty after applying --min_samples", file=sys.stderr)
        sys.exit(1)
    
    # Sort
    sort_cols = ['gene', 'position']
    if args.group_by:
        sort_cols = args.group_by + sort_cols
    entropy_df = entropy_df.sort_values(sort_cols)
    
    # Save
    print(f"\nSaving to: {args.output}")
    entropy_df.to_csv(args.output, index=False)
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Entropy records: {len(entropy_df):,}")
    print(f"Genes: {entropy_df['gene'].nunique()}")
    if args.group_by:
        for col in args.group_by:
            print(f"  {col} groups: {entropy_df[col].nunique()}")
    print(f"\nEntropy stats:")
    print(f"  Mean: {entropy_df['entropy'].mean():.3f} bits")
    print(f"  Median: {entropy_df['entropy'].median():.3f} bits")
    print(f"  Max: {entropy_df['entropy'].max():.3f} bits")
    print(f"  Min: {entropy_df['entropy'].min():.3f} bits")
    print(f"\nDone!")

if __name__ == "__main__":
    main()
