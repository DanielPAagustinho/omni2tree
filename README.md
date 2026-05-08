
# Welcome to Omni2Tree!!

Omni2Tree is an end-to-end phylogenomic workflow that builds maximum-likelihood phylogenetic trees directly from raw sequencing reads, without genome assembly or reliance on a single reference. It (i) constructs a reference database with OMA Standalone from coding sequences retrieved via NCBI accessions and/or extracted from local GTF/GFF3-annotated assemblies, (ii) processes read samples with optional deduplication and downsampling, and (iii) generates phylogenetic trees that place reference assemblies and read-derived consensus sequences on the same tree. For segmented viruses, Omni2Tree can be run on individual segments or selected subsets of segments, enabling users to investigate reassortment patterns or focus on specific genomic regions. The workflow also handles metagenomic samples containing multiple co-infecting pathogens, integrates metadata into an interactive HTML visualization through Omni2TreeView, computes per-position Shannon entropy from the merged MSAs with publication-ready plots, and includes a utility to bulk-download SRA runs or experiments.

## Table of Contents

- [Installation](docs/installation.md)
- [Quick start](#quick-start)
- [hRSV demo dataset](demo/README.md)
- [Running step 1: Creating the reference database](#running-step-1-creating-the-reference-database)
- [Running step 2: Processing sample reads and adding them to the read2tree folder](#running-step-2-processing-sample-reads-and-adding-them-to-the-read2tree-folder)
- [Running step 3: Getting the tree](#running-step-3-getting-the-tree)
- [Bulk Downloading read samples](#bulk-downloading-read-samples)
- [Shannon Entropy Analysis Pipeline for MSA Data](docs/entropy.md)

## Installation

See the full installation guide in [docs/installation.md](docs/installation.md).

After installing the dependencies, install Omni2Tree itself with:

```bash
git clone https://github.com/DanielPAagustinho/omni2tree.git
cd omni2tree
./install.sh /your/install/path
```

## Quick start

Minimal example (adjust paths to your data):

```bash
# Step 1: Create reference OMA database for r2t with NCBI and/or local references
o2t-step1 -i rsv_accessions.csv -g rsv_outgroups.csv -T 25 -o o2t_rsv &> rsv_long_step1.log

# Step 2: Map long nanopore reads to the reference using GNU Parallel
parallel -j 4 o2t-step2 -r {1} -o o2t_rsv -T 20 ::: \
  $(ls reads/*fastq* | sort) &>> "rsv_long_step2.log" &

# Step 2: Map short paired end reads to the reference
  parallel -j 4 o2t-step2 \
  -r {1} -t paired -map_op '"-ax sr"' -o o2t_rsv -T 20 ::: \
  $(ls reads/*_1.fastq* | sort) :::+ $(ls reads/*_2.fastq* | sort) &>> "rsv_short_step2.log" &

# Step 3: Create tree + visualization + entropy outputs (excluding reference sequences from the entropy analysis)
o2t-step3 -o o2t_rsv -m data/metadata.csv -l hRSV_main --seq_type aa --exclude_pattern "s0" -T 8
```
## hRSV demo dataset

A complete small hRSV demo is available in [demo/](demo/README.md).
It runs Step 1, Step 2, Step 3, metadata validation, entropy plots and
Omni2TreeView HTML generation.

Please, take into account that tree branch lengths/support values can vary between runs because MAFFT/IQ-TREE
use stochastic steps, so the demo should be validated with expected labels, row counts, generated files and log success patterns rather than exact Newick
identity.


## Running step 1: Creating the reference database

<details>
<summary>Click to expand/collapse</summary>

```bash
o2t-step1 -i rsv_accessions.csv -g rsv_outgroups.txt -T 25 -o virus2tree_rsv --temp_dir temp --debug &> def_rsv_long.log
```
To create the reference database, provide at least one reference source:

`-i`, `--input` (optional if `--local_assemblies` is provided): A file containing the NCBI accessions to be used for reference (`rsv_accessions.csv`).

`-L`, `--local_assemblies` (optional if `--input` is provided): A CSV manifest containing local assemblies and GTF/GFF3 annotations.

`-g`, `--outgroup` (optional but recommended): A file containing taxon(s) to be used as outgroups by OMA Standalone during orthologous group prediction (`rsv_outgroups`).  
If this file is not specified, OMA Standalone will use midpoint rooting, which is likely incorrect and will significantly affect hierarchical orthologous groups (HOGs) inferred by OMA.

### **Command Parameters**

| **Parameter**       | **Description** |
|--------------------|-----------------------------|
| `-i`, `--input`    | Optional CSV file with NCBI accessions. Required only when `--local_assemblies` is not provided. |
| `-L`, `--local_assemblies` | Optional CSV manifest with local assemblies to add to the same reference database as the NCBI accessions. Required only when `--input` is not provided. |
| `--local_features` | Comma-separated feature type(s) from local GTF/GFF3 column 3 to extract. **Default:** `CDS`. |
| `--local_group_by` | Optional single GTF/GFF3 attribute used to group selected local feature rows into one sequence unit. |
| `-g`, `--outgroup` | **Optional (recommended)** File with outgroup taxa used by OMA. |
| `-o`, `--o2t_out`       | Base output directory where all outputs are written. **Default:** current directory|
| `--temp_dir`       | Temporary directory (relative to `--o2t_out` or absolute). **Default:** `mktemp -d`|
| `-R`, `--resume`       | Skips taxa whose reference FASTA already exists in the `db` folder. If all taxa are already present, it resumes at Step 1.4. If the required cleaned files are already present and match `db/`, Step 1.4 is bypassed and the script resumes from the OMA Standalone run (Step 1.6). When run, it *removes* existing OMA output & read2tree directories. |
|`--og_min_fraction`| Keep only OGs present in at least this fraction of species (0–1). If omitted, all OGs are kept. |
| `-p, --use_mat_peptides`       | For NCBI accessions, download GBK files and use `mat_peptide` features instead of CDS features if at least one `mat_peptide` is found. |
| `-q, --use_mat_peptides_only`       | Same as `--use_mat_peptides`, except that if no `mat_peptide` feature is found, it does not download CDS features and skips that taxon. |
| `-T, --threads`   | Number of threads to use for OMA Standalone and the first step of read2tree. **Default:** `12`. |
| `--debug`         | Keeps temp directory with intermediate files. |
| `-h, --help`         | Show help. |

### **Accession File Format**

When `-i/--input` is used, the accession file must be a comma-separated values (CSV) text file, with the first line as the header. Each line represents a taxon/species/strain/label with associated accessions. The format varies depending on whether a five-letter code is included.

#### **Columns:**
1. **First column (required):** Taxon/species/label/strain name. Header: taxon (or taxa), species, label(s) or strain(s).
2. **Second column (optional):** Five-letter code. Must be exactly 5 alphanumeric characters. Header: code(s). If not provided, a random five-letter code for each taxon will be generated and saved in the file five_letter_taxon.tsv
3. **Third and onward (required):** One or more accession numbers (comma-separated) to obtain coding sequences. Accepts NCBI Nucleotide database accessions and assembly identifiers (GCF_/GCA_). Header: accession(s).

Commented lines starting with `#` are ignored.

For segmented viruses, include only the accessions corresponding to the segment(s) of interest if the goal is to analyze reassortment or focus on a subset of the genome.

#### **Example Input Files**

##### **With a five-letter code:**
```plaintext
STRAINS,CODE,accessions
influenza A virus California,INCFA,GCF_001343785.1
Influenza A Hong Kong, INHKA,GCF_000851145.1
ebola virus,EBOLA,GCA_034098425.1
Measles morbillivirus,MEAMO,GCF_000854845.1
Lyssavirus rabies,RABIE,GCF_000859625.1
Mammarenavirus lassaense,MAMMA,GCF_000851705.1
```

##### **Without a five-letter code:**
```plaintext
taxon,accession
chikungunya virus S-27,GCA_000854045.1
SARS-COV 2,NC_045512.2
Norovirus GI,NC_044853.1
Norovirus GIV,NC_044855.1
Norovirus GII,NC_044932.1
Norovirus GIII,NC_029645.1
Norovirus GV, NC_008311.1
```
> ⚠️ **Important**: Only the alphanumeric characters in the taxon column are considered for downstream processing. Taxon names and codes must be unique; duplicates are not allowed.

### **Outgroup File Format**

The outgroup file should contain taxa from the accession file to be used as outgroups. It can include one or more taxa.

#### **Example Outgroup File**

```plaintext
influenza A virus California
Influenza A Hong Kong
```

### **Local Assemblies Manifest**

`-L`, `--local_assemblies` adds homemade/local assemblies to the reference database built in step 1. Each row represents one reference taxon/label and must point to a DNA FASTA file plus a GTF/GFF3 annotation file. By default, Omni2Tree extracts features whose third GTF/GFF3 column is `CDS`, writes them to `db/{taxon}_cds_from_genomic.fna`, and then processes them together with the NCBI references if `-i/--input` is also provided.

If the main accession file has a code column, the local manifest must also have one:

```plaintext
taxon,code,dna,gff
Local RSV A,LRSA1,local_rsv_a.fasta,local_rsv_a.gff3
```

If the main accession file does not have a code column:

```plaintext
taxon,dna,gff
Local RSV A,local_rsv_a.fasta,local_rsv_a.gff3
```

If `-i/--input` is omitted, either local manifest format is accepted.

Use `--local_features` to select one or more feature types from GTF/GFF3 column 3, for example `--local_features CDS`, `--local_features mat_peptide`, or `--local_features CDS,mat_peptide`. By default, selected feature rows are not grouped: each selected row is treated as one sequence unit, even when multiple feature types are selected. Use `--local_group_by <attribute>` to group selected rows into final units by one attribute; a final group can contain rows from different selected feature types.

Each non-comment GTF/GFF3 row must have at least 9 tab-separated columns. Only include the feature rows you want to use. For example, if you want a subset of genes or segments, pre-filter the annotation so that only those rows remain. Local taxa must not duplicate taxa from the main accession file after Omni2Tree's alphanumeric taxon-name cleanup.

### **Output Files**

All outputs are placed within the `--o2t_out` directory

| **File**                      | **Description** |
|--------------------------------|------------------------------------------------------------------|
| `db/{taxon}_cds_from_genomic.fna` | Nucleotide FASTA files for each taxon with NCBI CDS/mature peptides or selected local GTF/GFF3 feature units. |
| `DB/{taxon}.fa`               | Amino acid FASTA files for each taxon, prepared for OMA Standalone and read2tree. |
| `dna_ref.fa`                  | Reference FASTA file with all cleaned nucleotide reference units from all taxa, prepared for read2tree. |
| `five_letter_taxon.tsv`        | Table linking taxa with five-letter codes. |
| `parameters.drw`              | Parameter file for the OMA run, modified according to the outgroup file and with the last 4 steps of OMA Standalone deactivated. |
| `Output/`                     | Folder containing the output from OMA Standalone. |
| `marker_genes/`               | Folder required by read2tree with all OMA orthologous group FASTA files, or only OGs passing `--og_min_fraction` when that option is used. |
| `stats/References_CDS/` | Directory containing per-reference sequence-unit count summaries generated in step 1. |
| `stats/References_CDS/cds_count_per_accession*` | Per-reference sequence counts and their distribution across NCBI and local references: `cds_count_per_accession.tsv`, `cds_count_per_accession_frequency.tsv` and `cds_count_per_accession_distribution.png`. |
| `stats/References_OGs/` | Directory containing reference OG summary outputs generated in step 1. |
| `stats/References_OGs/OG_genes.tsv` | Table with all features for each CDS from the OGs identified by OMA. |
| `stats/References_OGs/OG_genes-unique.tsv` | Summary table listing the OGs alongside its associated gene, protein, and the taxa in which it is found. |
| `stats/References_OGs/taxon_OG.tsv` | Table containing per-taxon summary: total CDS, missing protein_id, no-OG matches, and matched counts. |
| `stats/References_OGs/OG_taxa.tsv` | Summary of species coverage per OG and whether it is kept (only when `--og_min_fraction` is used). |
| `stats/References_OGs/og_taxa_per_og*` | Per-OG taxa counts and their distribution: `og_taxa_per_og.tsv`, `og_taxa_per_og_frequency.tsv` and `og_taxa_per_og_distribution.png`. |
| `stats/entropy/OG_genes_entropy.csv` | OG-to-gene mapping table generated from `stats/References_OGs/OG_genes.tsv` for the entropy workflow in step 3. |
| `O2T_RESULTS`     | Output of step 1 of read2tree. |

</details>

## Running step 2: Processing sample reads and adding them to the read2tree folder

<details>
<summary>Click to expand/collapse</summary>

After generating the reference database of orthologous groups, we proceed to add the sample reads.

```bash
#For long nanopore reads (Default for -t, --read_type is single and for --minimap2_options is "-ax map-ont")
parallel -j 4 o2t-step2 \
  -r {1} --dedup --downsample --coverage 250 --genome_size 15kb -o virus2tree_rsv -T 20 ::: \
  $(ls reads/*fastq* | sort) &>> "rsv_long_step2.log" &

#For paired end illumina reads
parallel -j 4 o2t-step2 \
  -r {1} {2} -t paired --minimap2_options '"-ax sr"' --dedup --downsample --coverage 250 --genome_size 15kb -o virus2tree_rsv -T 20 ::: \
  $(ls reads/*_1.fastq* | sort) :::+ $(ls reads/*_2.fastq* | sort) &>> "rsv_short_step2.log" &

```

### **Command Parameters**

| **Parameter**      | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `-r, --reads`     | **Required.** Input reads file(s) in `fastq` or `fastq.gz` format. If multiple files are provided and `--read_type` is not `paired`, they will be concatenated, assuming they belong to the same sample. |
| `-t, --read_type` | Generic read type: `single` or `paired`. If `paired`, two input files are required in `--reads`. **Default:** `single`.|
|`-map_op, --minimap2_options`| Options for minimap2 when mapping read set to the reference. Pass as a *single quoted string* (e.g., `--minimap2_options "-ax map-ont"`). Click [here](docs/recommended_presets.md) for suggested values. **Default:** `-ax map-ont`|
|`-o`, `--o2t_out`       | Base output directory that contains step 1 results and where step 2 writes outputs. **Default:** current directory.|
| `--temp_dir`      | Temporary directory (relative to `--o2t_out` or absolute). **Default:** `mktemp -d`. |
| `--stats_file`   | Name of the summary read statistics file. **Default:** `reads_statistics.tsv` | 
| `--dedup`        | Enables `czid-dedup` to remove duplicate reads. |
| `--dedup_l`      | Prefix length used for deduplication (requires `--dedup`). |
| `--downsample`   | Enables `rasusa` for read subsampling. It is required for all subsampling parameters |
| `--coverage`     | Minimum coverage for subsampling (integer or float: e.g., `250`, `0.1`). Requires `--genome_size`. |
| `--genome_size`  | Genome size for subsampling (integer or with a metric suffix: e.g., `15kb`, `4.1MB`). See [rasusa manual](https://github.com/mbhall88/rasusa?tab=readme-ov-file#genome-size) for more details. Requires `--coverage`. |
| `--num_bases`    | Target number of bases for subsampling (integer). Cannot to be used together with `--genome_size` and `--coverage`. |
| `--num_reads`    | Target number of reads for subsampling (integer). Cannot to be used together with `--genome_size` and `--coverage`. |
| `-T, --threads`   | Threads to use during step 2 of read2tree. **Default:** 4. |
| `--debug`        | Keeps temp directory with intermediate files. |
| `-h, --help`         | Show help. |

### **Output Files**

All outputs are placed within the `--o2t_out` directory

| **File**                       | **Description** |
|---------------------------------|------------------------------------------------------------------------------------------|
| `O2T_RESULTS`       | Directory used by steps 1 and 2 for read2tree output. |
| `reads_statistics.tsv`          | Summary of statistics for processed read samples, including initial state, deduplication, and downsampling. Reports the number of reads, average length, and total bases. |
| `temp/{sample}.fastq`           | Original reads. Uncompressed if initially compressed and concatenated if multiple input files were provided without `paired` option. |
| `temp/{sample}_dedup.fastq`     | Deduplicated reads. |
| `temp/{sample}_ds.fastq`        | Downsampled reads. |
| `temp/{sample}_dedup_ds.fastq`  | Deduplicated and downsampled reads. |

</details>

## Running step 3: Getting the tree

<details>
<summary>Click to expand/collapse</summary>

Finally, to obtain the phylogenetic tree and much more, run

```bash
o2t-step3 \
  -o virus2tree_rsv \
  -m data/metadata.csv \
  --seq_type aa \
  -T 8
```

This last step executes the following workflow:
1. Early metadata validation (`validate_metadata.py`)
2. `read2tree --step 3combine`
3. IQ-TREE inference
4. Tree/metadata preparation for visualization (`prepare_metadata_o2t_view.py`)
5. HTML view generation (`omni2treeview.py`)
6. Entropy table generation (`msa_to_position_table.py`, `calculate_entropy.py`)
7. Entropy plotting (`plot_entropy.R`)


### **Command Parameters**

| **Parameter** | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `-o`, `--o2t_out` | **Required.** Base output directory containing step1/step2 outputs (`O2T_RESULTS`, `marker_genes`, `dna_ref.fa`, `five_letter_taxon.tsv`). |
| `-m, --metadata` | **Required.** Metadata CSV used for validation, relabeling and entropy grouping/filtering. |
| `--seq_type` | Sequence mode: `aa` or `dna`. Selects the concatenated alignment used by IQ-TREE and entropy inputs. **Default:** `aa`. |
| `-T, --threads` | Threads for IQ-TREE. **Default:** `4`. |
| `--bootstrap` | IQ-TREE ultrafast bootstrap replicates. **Default:** `1000`. |
| `-r`, `--redo` | Allow rerunning step 3 in the same output directory. |
| `--temp_dir` | Optional temp directory. |
| `--debug` | Keep temporary files. |
| `-l, --label` | Optional visualization label. **Default:** `Omni2tree_Tree`. |
| `--reference_id` | Optional reference sequence ID or metadata label for entropy coordinates. When set, Step 3 maps alignment columns to ungapped reference positions and calculates entropy on those positions. |
| `--exclude_pattern` | Regex to exclude sample IDs during entropy table generation. Examples: `s0`, `RespiratorysyncytialvirusA`, `MinION_18_SRR33779449`, or a regex such as `^(MinION|NextSeq)_` to match multiple readsets. |
| `--filter_column` | Metadata column used to filter samples in entropy step 1 (requires `--filter_value`). |
| `--filter_value` | Metadata value kept in entropy step 1 (requires `--filter_column`). |
| `--group_by` | One or more metadata columns to group entropy calculations. |
| `--min_samples` | Minimum samples per position for entropy. |
| `--exclude_gaps` | Exclude gap character (`-`) before entropy calculation. |
| `--add_domain` | Optional domain CSV (`gene,domain,start,end`) to add domain ranges to per-gene entropy plots. Coordinates are interpreted in the entropy coordinate system: alignment positions by default, reference positions with `--reference_id`. The `gene` value must match the `gene` column in `stats/entropy/OG_genes_entropy.csv` generated in step 1 from `stats/References_OGs/OG_genes.tsv`. |
| `-h, --help` | Show help. |

### **Metadata Input Format (`-m`)**

`metadata.csv` must be a valid CSV with:
1. Header row
2. Type row (`character`, `date`, `numeric`, `integer`)
3. Data rows

Required columns (case-insensitive):
- `label`
- `accession`

Validation rules:
- `label` and `accession` cannot be empty.
- `label` and `accession` cannot contain commas.
- `label` must be unique.
- `accession` must be unique.
- Every reference label from `five_letter_taxon.tsv` must be represented in metadata (alphanumeric-clean matching).
- Every consensus readset detected in `O2T_RESULTS/*_all_cov.txt` must exist in metadata `label`.
- Do not add a column to mark rows as reference/readset; Step 3 detects this automatically.

Readset naming rule for `label`:
- Use the readset name exactly as generated by Step 2 (filename base without extension).
- For paired-end inputs, use the first mate readset name (e.g., `_1`/`R1` name).

### **Metadata Example**

```csv
label,accession,Genotype,Host
character,character,character,character
RespiratorysyncytialvirusA,NC_038235.1,Ref,Human
MinION_18_SRR33779449,SRR33779449,A,Human
NextSeq_22_SRR34709337_1,SRR34709337,A,Human
```

### **Main Outputs**

| **File/Directory** | **Description** |
|--------------------|-----------------|
| `O2T_RESULTS/tree_merge.nwk` | Newick tree generated by `read2tree --step 3combine`. |
| `O2T_RESULTS/concat_*_<aa\|dna>.phy.treefile` | IQ-TREE `treefile` generated from the concatenated alignment. |
| `O2T_RESULTS/concat_merge_view_<aa\|dna>.phy.treefile` | Tree rewritten by the wrapper for visualization after metadata relabeling. |
| `visualization/omni2treeview_<label>.html` | Interactive Omni2tree HTML output (`<label>` comes from `--label`). |
| `stats/entropy/<aa\|dna>_positions.csv` | Position-level table generated from the merged MSAs in the first entropy substep. |
| `stats/entropy/<aa\|dna>_entropy.csv` | Shannon entropy results. |
| `stats/entropy/entropy_all_genes.png` | Combined entropy plot across all genes. |
| `stats/entropy/entropy_<gene>.png` | Per-gene entropy plots; if `--add_domain` is provided, the matching domain ranges are highlighted here. |

</details>

## Bulk downloading read samples

<details>
<summary>Click to expand/collapse</summary>

To easily get read samples from SRA database, we developed the script `o2t-sra`. The purpose of this script is to facilitate the download and conversion to FASTQ of SRA IDs, whether they correspond to RUN (e.g., SRR, ERR, DRR) or EXPERIMENT (e.g., SRX, ERX, DRX).

Depending on the IDs type, the script proceeds differently:

* RUN mode: The script downloads each RUN individually and converts it to FASTQ.
* EXPERIMENT mode: The script first identifies which RUNs are associated with each EXPERIMENT, and then downloads each resulting RUN.

By default, the script retrieves metadata from the SRA database using `esearch` and `efetch`, then checks column 16 of runinfo to automatically determine whether each RUN is SINGLE or PAIRED. However, if you specify `-l` or `--layout` to force one layout (SINGLE or PAIRED), the script will skip the metadata-based check and apply the specified layout to all downloaded RUNs. However, if your inputs are EXPERIMENT IDs, the script will still need to retrieve runinfo in order to map each experiment to its associated RUNs, even if a forced layout is specified.

The output consists of FASTQ files corresponding to each RUN, renamed to include the species name and, in the case of EXPERIMENTs, also the experiment accession.

### **Command Parameters**

| **Parameter**      | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `-i, --input`     | **Required.** Input file containing SRA IDs with one taxon per line. The format should be `<species_name>,SRA_ID1,SRA_ID2,...`. |
| `-o, --outdir` | Directory where downloaded read files will be saved.	Default: current directory (`pwd`). |
| `-c, --chunk-size`   | Number of SRA IDs per chunk when fetching metadata using esearch and efetch. Default: `350`. |
| `-w, --sleep-secs`      | Number of seconds to sleep between chunked metadata requests to avoid overloading the NCBI server. Default: `1`. |
| `-l, --layout`   | Force sequencing layout (`SINGLE` or `PAIRED`) for all runs. When specified, skips metadata fetching. | 
|`-d, --debug`	| Keep per-species temporary directories and intermediate files. |   
| `-h, --help`        | Displays help information and exits. |

#### **Example Command**

```bash
o2t-sra -i sra_runs_rsv.csv --chunk-size 3 --sleep-secs 2 --outdir rsv_reads
```

It means to analyze the input file `sra_runs_rsv.csv`, extract the SRA IDs, fetch metadata in chunks of 3 SRA IDs at a time with a 2-second pause between chunks to avoid server overload, and download the resulting reads into the specified output directory `rsv_reads`.

> ⚠️ **Important**: You can adjust `--chunk-size` and `--sleep-secs` to avoid speed issues or overloading of the NCBI server

### **Input File Format**

The input file is a comma-separated values (CSV) text file. No header is required. Each line represents a taxon name (or any identifier of your choice) with one or more SRA IDs separated by commas. 

#### **Columns:**
1. **First column:** Species name (or any identifier). Spaces are allowed but should be avoided for simplicity. This is used for file naming and is sanitized to contain only alphanumeric characters
2. **Second column and onward:** One or more SRA IDs. All IDs in a line must be of the same type (RUN or EXPERIMENT). RUN: IDs starting with SRR, ERR, or DRR. EXPERIMENT: IDs starting with SRX, ERX, or DRX.

Commented lines starting with # are ignored.

#### **Example Input File**

```plaintext
SpeciesA,SRR123456,SRR123457
SpeciesB,SRR999999
# Example of a comment (this line is ignored)
SpeciesC,SRX000111
SpeciesD,ERX222333,ERX222334
```

### **Output Files**

At the end of processing, all generated fastq files are saved in the output directory:

For RUN mode (e.g., SRR123456 from Species A) it generates:
* `SpeciesA_SRR123456_1.fastq` and `SpeciesA_SRR123456_2.fastq` (if layout is PAIRED), or
* `SpeciesA_SRR123456.fastq` (if layout is SINGLE).

> **Note**: If the species name contains spaces (e.g., My Species), they will be deleted.

For EXPERIMENT mode (e.g., SRX000111 with RUN SRR000999 from Species A) it generates:
* `SpeciesA_SRX000111_SRR000999_1.fastq` and `SpeciesA_SRX000111_SRR000999_2.fastq` (if layout is PAIRED), or
* `SpeciesA_SRX000111_SRR000999.fastq` (if layout is SINGLE).

Also, a summary file with details about the SRA IDs downloaded is available:`{outdir}/summary_download.txt` (appends a per-species report + global totals).

At the end of execution, the script removes the directories containing the .sra and metadata files, leaving only the final FASTQ files and the summary file in the specified output directory.

</details>

## Citation

If you use Omni2Tree, please cite:

- **Omni2Tree**: Majidian, S., Chalco, A., Zheng, X., Webby, R. J., Bowman, A. S., Poulson, R. L., Nemeth, N. M., Sedlazeck, F. J., & Agustinho, D. P. (2026). "Rapid phylogenomic analysis for viral surveillance and metagenomic profiling with Omni2Tree" *bioRxiv*. https://doi.org/10.64898/2026.04.29.721707

---

## License

MIT License - Feel free to use and modify for your research.

---

## Support

For issues or questions:
1. Check this README thoroughly
2. Verify input file formats match examples
3. Test with a small subset of data first
4. Open an issue on GitHub with error messages and example data

---
