
from Bio import SeqIO
from Bio.Seq import Seq
import sys
from os import listdir
import os
from datetime import datetime

if sys.stdout.isatty():
    RED = "\033[1;31m"
    GREEN = "\033[1;32m"
    YELLOW = "\033[1;33m"
    NC = "\033[0m"
else:
    RED = GREEN = YELLOW = NC = ""

def log_info(message):
    print(f"{GREEN}[{datetime.now():%Y-%m-%d %H:%M:%S}] [INFO]{NC} {message}")

def log_warn(message):
    print(f"{YELLOW}[{datetime.now():%Y-%m-%d %H:%M:%S}] [WARN]{NC} {message}")

def log_error(message):
    print(f"{RED}[{datetime.now():%Y-%m-%d %H:%M:%S}] [ERROR]{NC} {message}", file=sys.stderr)

def read_fasta_files(input_folder):
    # for the same accession file generate the same codes
    files = sorted(listdir(input_folder))
    records_all = []
    file_names = [] 
    suffix = "_cds_from_genomic.fna"
    for file in files:
        if file.endswith(suffix):
            file_names.append(file)
            records = list(SeqIO.parse(input_folder + file, "fasta"))
            records_all.append(records)        
        else:
            log_info("Skipping file: " + str(input_folder + file) + " (not " + suffix + ")")
    if records_all:
        total_sequences = sum(len(r) for r in records_all)
        log_info(f"Found {len(file_names)} '{suffix}' files.")
        log_info(f"Total sequences: {total_sequences}")
        log_info(f"First file ({file_names[0]}) has {len(records_all[0])} sequences")
    else:
        log_info(f"No '{suffix}' files found in {input_folder}.")
    return file_names, records_all

def load_provided_codes(tsv_file):
    """Reads a TSV file with strain and 5-letter codes."""
    code_mapping = {}
    with open(tsv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                log_warn(f"Skipping invalid line in {tsv_file}: {line}")
                continue
            strain, code = parts
            code_mapping[strain] = code
    log_info(f"Loaded {len(code_mapping)} strain-code pairs from {tsv_file}.")
    return code_mapping

def create_five_letter(file_names):
    """Generates 5-letter codes for each file."""
    five_letter_species_dic = {}
    for idx, file_name in enumerate(file_names):
        five_letter_species = "s" + str(idx).zfill(4) 
        five_letter_species_dic[file_name.replace("_cds_from_genomic.fna","")] = five_letter_species
    log_info("Generated automatic 5-letter codes.")
    return five_letter_species_dic

def clean_translate(records, species_fivelet, strain):
    """Cleans and translates nucleotide sequences."""
    records_nuc = []
    records_aa = []
    for record in records:
        sequence = record.seq
        remainder = len(sequence) % 3
        if remainder != 0:
            sequence += Seq('N' * (3 - remainder))
            record.seq = sequence
        #Delete trouble symbols
        id_old = str(record.id).replace("_", "").replace(".", "")
        strain =''.join(c for c in strain if c.isalnum()) #really necessary? this has already been solved by fetch_data
        id_new = species_fivelet+"||"+strain + "||" + id_old
        
        nuc_seq = SeqIO.SeqRecord(sequence, id=id_new, description="cleaned for r2t", name=id_new)
        protein_seq = SeqIO.SeqRecord(sequence.translate(), id=id_new, description="cleaned for r2t", name=id_new)
        
        records_nuc.append(nuc_seq)
        records_aa.append(protein_seq)
    
    log_info(f"Cleaned nucleotide and protein sequences for {species_fivelet}.")
    return records_nuc, records_aa

if __name__ == '__main__':
    if len(sys.argv) < 2:
        log_info("Usage: python script.py <input_folder_fna> [resuming_download] [strain_code_tsv]")
        sys.exit(1)

    input_folder_fna = sys.argv[1].rstrip('/') + "/"
    resuming_download = sys.argv[2].lower() == "true" if len(sys.argv) > 2 else False    # Read input FASTA files
    strain_code_tsv = sys.argv[3] if len(sys.argv) > 3 else None
    
    file_names, records_all = read_fasta_files(input_folder_fna)

    # Determine whether to use provided codes or generate them
    if strain_code_tsv and os.path.exists(strain_code_tsv):
        five_letter_species_dic = load_provided_codes(strain_code_tsv)
        # Ensure all files have a corresponding code
        for file_name in file_names:
            strain = file_name.replace("_cds_from_genomic.fna", "")
            if strain not in five_letter_species_dic:
                log_error(f"No 5-letter code provided for strain '{strain}'.")
                sys.exit(1)
    else:
        five_letter_species_dic = create_five_letter(file_names)
        # Output clean_five_letter_species.tsv only if no TSV was provided
        output_five_letter_tsv = "five_letter_taxon.tsv"
        with open(output_five_letter_tsv, "w") as file_out:
            for species_name, five_letter in five_letter_species_dic.items():
                file_out.write(f"{species_name}\t{five_letter}\n")
        log_info(f"Generated codes written to {output_five_letter_tsv}.")

    # Create output folder for amino acid sequences
    folder_aa = "DB"
    if not os.path.exists(folder_aa):
        os.makedirs(folder_aa)
    elif resuming_download:
        log_warn(f"The folder '{folder_aa}' already exists. New sequence(s) will be added to it.")
    else:
        log_error(f"The folder '{folder_aa}' already exists. Please remove it.")
        sys.exit(1)

    # Process files and clean/translate sequences
    records_nuc_all_clean = []
    for idx in range(len(file_names)):
        file_name = file_names[idx]
        records = records_all[idx]
        strain = file_name.replace("_cds_from_genomic.fna", "")
        species_fivelet = five_letter_species_dic[strain]
        
        records_nuc, records_aa = clean_translate(records, species_fivelet,strain)
        SeqIO.write(records_aa, f"{folder_aa}/{strain}.fa", "fasta")
        records_nuc_all_clean += records_nuc

    # Write all nucleotide sequences to dna_ref.fa
    SeqIO.write(records_nuc_all_clean, "dna_ref.fa", "fasta")
    log_info(f"Processed {len(file_names)} files. Cleaned nucleotide sequences saved to 'dna_ref.fa'.")
    log_info(f"Amino acid sequences saved to '{folder_aa}'.")
