#!/usr/bin/env python3
"""
Convert MSA files to a position table for entropy analysis.
Handles both amino acid and nucleotide sequences in FASTA/Phylip format.

Usage:
    python msa_to_position_table.py --msa_dir MSA/hepc/AA \
                                     --og_table hepc_ogs.csv \
                                     --output hepc_positions.csv \
                                     --metadata metadata.csv \
                                     --seq_type AA

Author: Daniel Agustinho
"""

import argparse
import pandas as pd
from pathlib import Path
from Bio import AlignIO
import sys
import re


METADATA_TYPE_VALUES = {"character", "numeric", "integer", "date", "factor", "logical"}

def read_phylip_msa(filepath):
    """
    Read MSA in Phylip format (Read2Tree output format).
    Handles the specific format with sequence length header.
    """
    try:
        alignment = AlignIO.read(filepath, "phylip-relaxed")
        return alignment
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return None

def parse_msa_to_dataframe(msa_file, og_name, gene_name, seq_type):
    """
    Parse MSA file and convert to long-format dataframe.
    
    Parameters:
    -----------
    msa_file : Path
        Path to MSA file
    og_name : str
        Ortholog group identifier (e.g., OG1)
    gene_name : str
        Gene/peptide name from mapping table
    seq_type : str
        'AA' for amino acids or 'DNA' for nucleotides
    
    Returns:
    --------
    pd.DataFrame with columns: sample_id, position, character, og, gene, seq_type
    """
    alignment = read_phylip_msa(msa_file)
    
    if alignment is None:
        return pd.DataFrame()
    
    records = []
    
    for record in alignment:
        sample_id = record.id.strip()
        sequence = str(record.seq)
        
        for position, char in enumerate(sequence, start=1):
            records.append({
                'sample_id': sample_id,
                'position': position,
                'character': char,
                'og': og_name,
                'gene': gene_name,
                'seq_type': seq_type
            })
    
    return pd.DataFrame(records)

def load_og_mapping(og_table_file):
    """
    Load OG to gene name mapping.
    Expected format: OG,peptide (or OG,gene)
    """
    og_table = pd.read_csv(og_table_file)
    
    # Handle different possible column names
    if 'peptide' in og_table.columns:
        og_table = og_table.rename(columns={'peptide': 'gene'})
    
    return dict(zip(og_table['OG'], og_table['gene']))


def clean_alnum(text):
    """Match step1 cleaning behavior: keep only alphanumeric characters."""
    return re.sub(r'[^A-Za-z0-9]', '', '' if pd.isna(text) else str(text))


def normalize_msa_id(sample_id):
    """Remove readset suffixes added by mapping/merging stages."""
    normalized = str(sample_id).strip()
    return re.sub(r'(_R[12]|_[12])$', '', normalized)


def maybe_drop_metadata_type_row(metadata):
    """
    Drop type row (character/numeric/...) if present as first data row.
    """
    if metadata.empty:
        return metadata

    first_row = metadata.iloc[0].astype(str).str.strip().str.lower()
    if first_row.apply(lambda x: x in METADATA_TYPE_VALUES).all():
        print("Detected metadata type row, excluding it from matching/filtering")
        return metadata.iloc[1:].reset_index(drop=True)
    return metadata


def load_five_letter_codes(five_letter_file):
    """
    Load five_letter_taxon.tsv and index by cleaned taxon label.
    Format expected: taxon<tab>code
    """
    labelkey_to_code = {}
    with open(five_letter_file, 'r', encoding='utf-8') as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                print(
                    f"Warning: Invalid five_letter line {line_no}, expected 2 tab-separated columns: {line}",
                    file=sys.stderr,
                )
                continue
            taxon, code = parts[0].strip(), parts[1].strip()
            label_key = clean_alnum(taxon)
            if not label_key:
                continue
            labelkey_to_code[label_key] = code
    return labelkey_to_code


def build_metadata_label_lookup(metadata, match_column, five_letter_map):
    """
    Build mapping from normalized MSA sample IDs to output labels.
    Matching is label-first (default) and additionally uses five-letter codes
    derived from cleaned labels for references.
    """
    if match_column not in metadata.columns:
        raise ValueError(f"Column '{match_column}' not found in metadata")

    output_label_column = 'label' if 'label' in metadata.columns else match_column
    lookup = {}
    collisions = 0

    for _, row in metadata.iterrows():
        match_value = row[match_column]
        if pd.isna(match_value):
            continue
        match_value = str(match_value).strip()
        if not match_value:
            continue

        output_label = row[output_label_column]
        if pd.isna(output_label):
            continue
        output_label = str(output_label).strip()
        if not output_label:
            continue

        candidate_ids = set()

        # Direct label/readset match (after suffix normalization)
        candidate_ids.add(normalize_msa_id(match_value))

        cleaned_label = clean_alnum(match_value)
        if cleaned_label:
            candidate_ids.add(normalize_msa_id(cleaned_label))

        # Reference match through five-letter code (taxon label -> sXXXX)
        code = five_letter_map.get(cleaned_label)
        if code:
            candidate_ids.add(normalize_msa_id(code))

        for candidate in candidate_ids:
            if not candidate:
                continue
            existing = lookup.get(candidate)
            if existing is not None and existing != output_label:
                collisions += 1
                continue
            lookup[candidate] = output_label

    if collisions > 0:
        print(
            f"Warning: {collisions} metadata collisions detected while building ID->label map; kept first occurrence",
            file=sys.stderr,
        )

    return lookup

def main():
    parser = argparse.ArgumentParser(
        description='Convert MSA files to position table for entropy analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process amino acid MSAs
  python msa_to_position_table.py --msa_dir MSA/hepc/AA --og_table hepc_ogs.csv \\
                                   --output hepc_aa_positions.csv --seq_type AA

  # Process DNA MSAs with metadata filtering
  python msa_to_position_table.py --msa_dir MSA/hepc/DNA --og_table hepc_ogs.csv \\
                                   --output hepc_dna_positions.csv --seq_type DNA \\
                                   --metadata metadata.csv --filter_column genotype \\
                                   --filter_value GT1
        """
    )
    
    parser.add_argument('--msa_dir', required=True, type=Path,
                        help='Directory containing MSA files (e.g., MSA/hepc/AA)')
    parser.add_argument('--og_table', required=True, type=Path,
                        help='CSV file mapping OG to gene/peptide names')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output CSV file for position table')
    parser.add_argument('--seq_type', required=True, choices=['AA', 'DNA'],
                        help='Sequence type: AA (amino acid) or DNA (nucleotide)')
    parser.add_argument('--metadata', type=Path, default=None,
                        help='Optional metadata CSV used for filtering and label mapping')
    parser.add_argument('--metadata_match_column', type=str, default='label',
                        help='Metadata column used to match MSA samples (default: label)')
    parser.add_argument('--five_letter', type=Path, default=None,
                        help='five_letter_taxon.tsv used to map cleaned metadata labels to 5-letter sample IDs')
    parser.add_argument('--filter_column', type=str, default=None,
                        help='Column name in metadata to filter by')
    parser.add_argument('--filter_value', type=str, default=None,
                        help='Value to filter for in filter_column')
    parser.add_argument('--exclude_pattern', type=str, default=None,
                        help='Pattern to exclude from sample IDs (e.g., s0 for reference strains). If not specified, all samples are included.')
    parser.add_argument('--include_all', action='store_true',
                        help='Include all samples (override exclude_pattern)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.msa_dir.exists():
        print(f"Error: MSA directory not found: {args.msa_dir}", file=sys.stderr)
        sys.exit(1)
    
    if not args.og_table.exists():
        print(f"Error: OG table not found: {args.og_table}", file=sys.stderr)
        sys.exit(1)
    
    # Load OG mapping
    print(f"Loading OG mapping from {args.og_table}...")
    og_mapping = load_og_mapping(args.og_table)
    print(f"Found {len(og_mapping)} OG entries")
    
    # Load metadata if provided
    metadata = None
    metadata_label_lookup = None
    if args.metadata:
        if not args.metadata.exists():
            print(f"Error: Metadata file not found: {args.metadata}", file=sys.stderr)
            sys.exit(1)
        print(f"Loading metadata from {args.metadata}...")
        metadata = pd.read_csv(args.metadata)
        metadata = maybe_drop_metadata_type_row(metadata)
        
        # Apply filtering if specified
        if args.filter_column and args.filter_value:
            if args.filter_column not in metadata.columns:
                print(f"Error: Column {args.filter_column} not found in metadata", file=sys.stderr)
                sys.exit(1)
            metadata = metadata[metadata[args.filter_column] == args.filter_value]
            print(f"Filtered to {len(metadata)} samples where {args.filter_column}={args.filter_value}")

        # Build metadata lookup with label-based matching by default
        five_letter_map = {}
        if args.metadata_match_column == 'label':
            if args.five_letter is None:
                print(
                    "Error: --five_letter is required when matching metadata by label",
                    file=sys.stderr,
                )
                sys.exit(1)
            if not args.five_letter.exists():
                print(f"Error: five_letter file not found: {args.five_letter}", file=sys.stderr)
                sys.exit(1)
            print(f"Loading five-letter mapping from {args.five_letter}...")
            five_letter_map = load_five_letter_codes(args.five_letter)
            print(f"Loaded {len(five_letter_map)} cleaned label -> 5-letter code mappings")

        try:
            metadata_label_lookup = build_metadata_label_lookup(
                metadata=metadata,
                match_column=args.metadata_match_column,
                five_letter_map=five_letter_map,
            )
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

        print(f"Prepared {len(metadata_label_lookup)} sample-id mappings from metadata")
    
    # Process all MSA files
    all_positions = []
    
    msa_files = sorted(args.msa_dir.glob("OG*.fa"))
    
    if not msa_files:
        print(f"Error: No OG*.fa files found in {args.msa_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nProcessing {len(msa_files)} MSA files...")
    
    for msa_file in msa_files:
        og_name = msa_file.stem  # e.g., OG1
        
        if og_name not in og_mapping:
            print(f"Warning: {og_name} not found in OG mapping table, skipping...", file=sys.stderr)
            continue
        
        gene_name = og_mapping[og_name]
        print(f"  Processing {og_name} ({gene_name})...")
        
        df = parse_msa_to_dataframe(msa_file, og_name, gene_name, args.seq_type)
        
        if df.empty:
            print(f"    Warning: No data extracted from {msa_file}", file=sys.stderr)
            continue
        
        # Filter by exclude pattern
        if args.exclude_pattern:
            before_count = df['sample_id'].nunique()
            df = df[~df['sample_id'].str.contains(args.exclude_pattern)]
            after_count = df['sample_id'].nunique()
            print(f"    Excluded {before_count - after_count} samples matching '{args.exclude_pattern}'")
        
        # Clean sample IDs (remove _R1/_R2 and _1/_2 suffixes)
        df['sample_id'] = df['sample_id'].apply(normalize_msa_id)
        
        # Filter by metadata if provided
        if metadata is not None:
            before_count = df['sample_id'].nunique()
            df = df[df['sample_id'].isin(metadata_label_lookup)]
            after_count = df['sample_id'].nunique()
            if not df.empty:
                # Output keeps human-readable labels exactly as provided in metadata.
                df['sample_id'] = df['sample_id'].map(metadata_label_lookup)
            print(f"    Filtered to {after_count}/{before_count} samples using metadata")
        
        all_positions.append(df)
        print(f"    Extracted {len(df)} position records from {df['sample_id'].nunique()} samples")
    
    if not all_positions:
        print("Error: No data was processed", file=sys.stderr)
        sys.exit(1)
    
    # Combine all dataframes
    print("\nCombining all position tables...")
    final_df = pd.concat(all_positions, ignore_index=True)
    final_df = final_df.rename(columns={'sample_id': 'label'})
    
    # Save output
    print(f"Saving to {args.output}...")
    final_df.to_csv(args.output, index=False)
    
    print(f"\n=== Summary ===")
    print(f"Total records: {len(final_df):,}")
    print(f"Unique samples: {final_df['label'].nunique()}")
    print(f"Genes processed: {final_df['gene'].nunique()}")
    print(f"Sequence type: {args.seq_type}")
    print(f"\nOutput saved to: {args.output}")

if __name__ == "__main__":
    main()
