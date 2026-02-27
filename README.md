
# Welcome to Omni2tree!!

This new version of read2tree enables the creation of a reference database via OMA Standalone using as input coding sequences from NCBI assemblies. The final tree combines the presence of assemblies (that is, the reference) and the read samples. It also supports read deduplication with czid-dedup and downsampling with rasusa, among other functionalities.

## Table of Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Running step 1: Creating the reference database](#running-step-1-creating-the-reference-database)
- [Running step 2: Processing sample reads and adding them to the read2tree folder](#running-step-2-processing-sample-reads-and-adding-them-to-the-read2tree-folder)
- [Running step 3: Getting the tree](#running-step-3-getting-the-tree)
- [Downloading read samples](#downloading-read-samples)
- [Shannon Entropy Analysis Pipeline for MSA Data](#shannon-entropy-analysis-pipeline-for-msa-data)

## Installation

### Dependencies

This software relies on four external tools: [OMA Standalone](https://omabrowser.org/standalone/), [Rasusa](https://github.com/mbhall88/rasusa?tab=readme-ov-file#install), [czid-dedup](https://github.com/chanzuckerberg/czid-dedup?tab=readme-ov-file#installation), and [Read2Tree](https://github.com/DessimozLab/read2tree/tree/minimap2?tab=readme-ov-file#installation). It assumes all programs are in your Conda environment or `PATH`. 
Below are two general ways to install all the required dependencies. For more details, please visit the respective web pages.

### 1. Installation with Conda (not possible for all dependencies)

[Conda](https://docs.anaconda.com/miniconda/) is a package manager that allows you to install all dependencies quickly and easily.

```bash
conda create -n my_env python=3.10.8 -y 
conda activate my_env 
conda install -c bioconda rasusa read2tree sra-tools entrez-direct -y
```

**Notes:** 
* OMA standalone and czid-dedup are not available via Conda. Please, follow the "Installation from source" instructions below.
* The Conda version of read2tree does not include the metagenomics branch. If you need this branch, follow the "Installation from source" instructions.

### 2. Installation from source

If you prefer to install the tools manually from their source code, use the following commands:

**OMA Standalone**

```bash
## Download the last version, in this example is 2.6.0
wget -O oma.tgz https://omabrowser.org/standalone/OMA.2.6.0.tgz 
tar xvzf oma.tgz 
cd OMA.2.6.0

## Below choose your install path, if not OMA will be installed in /usr/local/OMA (you might need to use sudo in this case)
./install.sh /your/install/path

## This step is IMPORTANT:
## After installation, make sure the bin folder of OMA is in your PATH variable. For that, edit your shell configuration file (`~/.bashrc`, `~/.zshrc`, etc.)
echo 'export PATH=$PATH:/your/install/path/OMA/bin' >> ~/.bashrc 
source ~/.bashrc
```

**Rasusa**

```bash
## When rasusa is downloaded it is automatically added to your PATH
curl -sSL rasusa.mbh.sh | sh
```

**czid-dedup**

Take into account that `czid-dedup` requires [rust/cargo](https://www.rust-lang.org/tools/install) for compilation

```bash
git clone https://github.com/chanzuckerberg/czid-dedup.git 
cd czid-dedup 
cargo build --release 

#Make sure that the release directory is in your PATH variable
echo 'export PATH=$PATH:your/install/path/czid-dedup/target/release' >> ~/.bashrc 
source ~/.bashrc
```

**Read2Tree**

```bash
## Create conda env
conda create -n r2t python=3.10.8 -y 
conda activate r2t

## Get required python packages
conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml filelock -y
conda install -c bioconda dendropy pysam -y

## Install required softwares
conda install -c bioconda mafft iqtree minimap2 samtools -y

## Clone read2tree
git clone https://github.com/DessimozLab/read2tree.git 
cd read2tree 
python setup.py install
## read2tree will be placed in the default bin folder of your Conda installation
```

**SRA Toolkit**

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz 
tar -xvzf sratoolkit.current-ubuntu64.tar.gz

## Add executable to your path (using your own version, in this case is 3.2.0)
echo 'export PATH="$PATH:/your/install/path/sratoolkit.3.2.0-ubuntu64/bin"' >> ~/.bashrc 
source ~/.bashrc
```

**Entrez direct**

```bash
## Get the scripts and download them in an "edirect folder" in the user's home directory
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
source ~/.bashrc
```

### 3. Verify Installation

To verify that all tools are correctly installed and available in your Conda environment or `PATH`, run the following command:

```bash
command -v oma && command -v rasusa && command -v czid-dedup && command -v read2tree && command -v fasterq-dump && command -v esearch
```

### Omni2tree installation

You can set up Omni2tree cloning the repo and running the installer:

```bash
git clone https://github.com/DanielPAagustinho/virus2tree.git
cd virus2tree
./install.sh /your/install/path
```

The installation script creates symlinks to the shell entry points.
* If you omit the install path, symlinks go to /usr/local/bin (may require sudo).
* If you use a custom path, make sure it’s in your PATH.

Finally, check your installation with:

```bash
which o2t-step1 && o2t-step1 --help
which o2t-step2 && o2t-step2 --help
which o2t-sra   && o2t-sra   --help
```

## Quick start

Minimal example (adjust paths to your data):

```bash
# Step 1: Create reference OMA database for r2t with NCBI accessions from RSV
o2t-step1 -i rsv_accessions.csv -g rsv_outgroups.txt -T 25 --o2t_out virus2tree_rsv &> rsv_long_step1.log

# Step 2: Map long nanopore RSV reads to the reference using GNU Parallel
parallel -j 4 o2t-step2 -r {1} --o2t_out virus2tree_rsv -T 20 ::: \
  $(ls reads/*fastq* | sort) &>> "rsv_long_step2.log" &

# Step 2: Map short paired end RSV reads to the reference
  parallel -j 4 o2t-step2 \
  -r {1} -t paired -map_op '"-ax sr"' --o2t_out virus2tree_rsv -T 20 ::: \
  $(ls reads/*_1.fastq* | sort) :::+ $(ls reads/*_2.fastq* | sort) &>> "rsv_short_step2.log" &

# Step 3: Create tree + visualization + entropy outputs
o2t-step3 --o2t_out virus2tree_rsv -m data/metadata.csv --seq_type aa -T 8
```

## Running step 1: Creating the reference database

<details>
<summary>Click to expand/collapse</summary>

```bash
o2t-step1 -i rsv_accessions.csv -g rsv_outgroups.txt -T 25 --o2t_out virus2tree_rsv --temp_dir temp --debug &> def_rsv_long.log
```
To create the reference database, two key input files are required:

`-i`, `--input` (Required): A file containing the NCBI accessions to be used for reference (`rsv_accessions.csv`).

`-g`, `--outgroup` (optional but recommended): A file containing taxon(s) to be used as outgroups by OMA Standalone during orthologous group prediction (`rsv_outgroups`).  
If this file is not specified, OMA Standalone will use midpoint rooting, which is likely incorrect and will significantly affect hierarchical orthologous groups (HOGs) inferred by OMA.

### **Command Parameters**

| **Parameter**       | **Description** |
|--------------------|-----------------------------|
| `-i`, `--input`    | **Required.** CSV file with NCBI accessions. |
| `-g`, `--outgroup` | **Optional (recommended)** File with outgroup taxa used by OMA. |
| `--o2t_out`       | Base output directory where all outputs are written. **Default:** current directory|
| `--temp_dir`       | Temporary directory (relative to `--o2t_out` or absolute). **Default:** `mktemp -d`|
| `--resume`       |Skips taxa already downloaded from NCBI into the `db` folder. If all taxa were already downloaded, it resumes at Step 1.4. Moreover, if the required files are already present, Step 1.4 is bypassed and the script practically resumes from the OMA Standalone run (Step 1.6). When run, it *removes* existing OMA output & read2tree directories. |
|`--og_min_fraction`| Keep only OGs present in at least this fraction of species (0–1). If omitted, all OGs are kept. |
| `-p, --use_mat_peptides`       | Download GBK files for each taxon's accession(s) and uses the mat_peptide features instead of CDS features if at least one mat_peptide is found. |
| `-q, --use_mat_peptides_only`       | Same as --use_mat_peptides, except that if no mat_peptide feature is found, it does not download CDS features and simply skips that taxon. |
| `-T, --threads`   | Number of threads to use for Oma Standalone and the first step of read2tree. |
| `--debug`         | Keeps temp directory with intermediate files. |
| `-h, --help`         | Show help. |

### **Accession File Format**

The accession file must be a comma-separated values (CSV) text file, with the first line as the header. Each line represents a taxon/species/strain with associated accessions. The format varies depending on whether a five-letter code is included.

#### **Columns:**
1. **First column (required):** Taxon/species/label/strain name. Header: taxon (or taxa), species, label(s) or strain(s).
2. **Second column (optional):** Five-letter code. Must be exactly 5 alphanumeric characters. Header: code(s). If not provided, a random five-letter code for each taxon will be generated and saved in the file five_letter_taxon.tsv
3. **Third and onward (required):** One or more accession numbers (comma-separated) to obtain coding sequences. Accepts NCBI Nucleotide database accessions and assembly identifiers (GCF_/GCA_). Header: accession(s).

Commented lines starting with `#` are ignored.

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

### **Output Files**

All outputs are placed within the `--o2t_out` directory

| **File**                      | **Description** |
|--------------------------------|------------------------------------------------------------------|
| `db/{taxon}_cds_from_genomic.fna` | Nucleotide FASTA files for each taxon with CDS (or mature peptides if selected) retrieved from NCBI. |
| `DB/{taxon}.fa`               | Amino acid FASTA files for each taxon, prepared for running with OMA Standalone and read2tree. |
| `dna_ref.fa`                  | Reference FASTA file with all nucleotide CDS from all taxa, prepared to be used as input for read2tree. |
| `five_letter_taxon.tsv`        | Table linking taxa with five-letter codes. |
| `parameters.drw`              | Parameter file for the OMA run, modified according to the outgroup file and with the last 4 steps of OMA Standalone deactivated. |
| `Output/`                     | Folder containing the output from OMA Standalone. |
| `marker_genes/`               | Folder required by read2tree with the orthologous groups (OGs) generated by OMA Standalone (contents of `Output/OrthologousGroupsFasta`). |
| `stats/cds_count_per_accession*`                | Per-assembly CDS counts and their distribution across all downloaded assemblies: `cds_count_per_accession.tsv`, `cds_count_per_accession_frequency.tsv` and `cds_count_per_accession.png` |
| `stats/OG_genes.tsv`                | Table with all features for each CDS from the OGs identified by OMA. |
| `stats/OG_genes-unique.tsv`         | Summary table listing the OGs alongside its associated gene, protein, and the taxa in which it is found. |
|`stats/taxon_OG.tsv`| Table containing per-taxon summary: total CDS, missing protein_id, no-OG matches, and matched counts. |
|`stats/OG_taxa.tsv	`|Summary of species coverage per OG and whether it is kept (only when --og_min_fraction is used)|
| `O2T_RESULTS`     | Output of step 1 of read2tree. |

</details>

## Running step 2: Processing sample reads and adding them to the read2tree folder

<details>
<summary>Click to expand/collapse</summary>

After generating the reference database of orthologous groups, we proceed to add the sample reads.

```bash
#For long nanopore reads (Default for -t, --read_type is single and for --minimap2_options is "-ax mp-ont")
parallel -j 4 o2t-step2 \
  -r {1} --dedup --downsample --coverage 250 --genome_size 15kb --o2t_out virus2tree_rsv -T 20 ::: \
  $(ls reads/*fastq* | sort) &>> "rsv_long_step2.log" &

#For paired end illumina reads
parallel -j 4 o2t-step2 \
  -r {1} {2} -t paired --minimap2_options '"-ax sr"' --dedup --downsample --coverage 250 --genome_size 15kb --o2t_out virus2tree_rsv -T 20 ::: \
  $(ls reads/*_1.fastq* | sort) :::+ $(ls reads/*_2.fastq* | sort) &>> "rsv_short_step2.log" &

```

### **Command Parameters**

| **Parameter**      | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `-r, --reads`     | **Required.** Input reads file(s) in `fastq` or `fastq.gz` format. If multiple files are provided and `--read_type` is not `paired`, they will be concatenated, assuming they belong to the same sample. |
| `-t, --read_type` | Generic read type: `single` or `paired`. If `paired`, two input files are required in `--reads`. **Default:** `single`.|
|`-map_op, --minimap2_options`| Options for minimap2 when mapping read set to the reference. Pass as a *single quoted string* (e.g., `--minimap2_options "-ax map-ont"`). Click [here](docs/recommended_presets.md) for suggested values. **Default:** `-ax map-ont`|
|`--o2t_out`       | Base output directory that contains step 1 results and where step 2 writes outputs. **Default:** current directory.|
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

Run the Omni2tree wrapper `o2t-step3` (not raw `read2tree` directly).  
This wrapper executes the full Step 3 workflow:
1. Early metadata validation (`validate_metadata.py`)
2. `read2tree --step 3combine`
3. IQ-TREE inference
4. Tree/metadata preparation for visualization (`prepare_metadata_o2t_view.py`)
5. HTML view generation (`omni2treeview.py`)
6. Entropy table generation (`msa_to_position_table.py`, `calculate_entropy.py`)
7. Entropy plotting (`plot_entropy.R`)

### **Example Command**

```bash
o2t-step3 \
  --o2t_out virus2tree_rsv \
  -m data/metadata.csv \
  --seq_type aa \
  -T 8
```

### **Command Parameters**

| **Parameter** | **Description** |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------|
| `--o2t_out` | **Required.** Base output directory containing step1/step2 outputs (`O2T_RESULTS`, `marker_genes`, `dna_ref.fa`, `five_letter_taxon.tsv`). |
| `-m, --metadata` | **Required.** Metadata CSV used for validation, relabeling and entropy grouping/filtering. |
| `--seq_type` | Sequence mode: `aa` or `dna`. Selects the concatenated alignment used by IQ-TREE and entropy inputs. **Default:** `aa`. |
| `-T, --threads` | Threads for IQ-TREE. **Default:** `4`. |
| `--bootstrap` | IQ-TREE ultrafast bootstrap replicates. **Default:** `1000`. |
| `--temp_dir` | Optional temp directory. |
| `--debug` | Keep temporary files. |
| `-l, --label` | Optional visualization label. **Default:** `Omni2tree_Tree`. |
| `--exclude_pattern` | Regex to exclude sample IDs during entropy table generation. |
| `--filter_column` | Metadata column used to filter samples in entropy step 1 (requires `--filter_value`). |
| `--filter_value` | Metadata value kept in entropy step 1 (requires `--filter_column`). |
| `--group_by` | One or more metadata columns to group entropy calculations. |
| `--min_samples` | Minimum samples per position for entropy. |
| `--exclude_gaps` | Exclude gap character (`-`) before entropy calculation. |
| `--add_domain` | Optional domain CSV (`gene,domain,start,end`) validated before plotting (currently validation-only; not rendered by `plot_entropy.R`). |
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
| `O2T_RESULTS/tree_merge.nwk` | Read2Tree merged tree from `3combine` (used as fallback if IQ-TREE `.treefile` is unavailable). |
| `O2T_RESULTS/output/concat_*_<aa|dna>.phy.treefile` or `O2T_RESULTS/concat_*_<aa|dna>.phy.treefile` | IQ-TREE inferred tree. |
| `O2T_RESULTS/concat_merge_view_<aa|dna>.phy.treefile` | Relabeled tree used for visualization. |
| `visualization/metadata_o2t_view.csv` | Viewer-ready metadata table. |
| `visualization/omni2treeview_<label>.html` | Interactive Omni2tree HTML output (`<label>` comes from `--label`, default `Omni2tree_Tree`). |
| `stats/entropy/<aa|dna>_positions.csv` | Position-level table from MSA. |
| `stats/entropy/<aa|dna>_entropy.csv` | Shannon entropy results. |
| `stats/entropy/*.png` | Entropy plots (combined + per gene). |

</details>

## Downloading read samples

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
./o2t-sra -i sra_runs_rsv.csv --chunk-size 3 --sleep-secs 2 --outdir rsv_reads
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

## Shannon Entropy Analysis Pipeline for MSA Data
<details>
<summary><strong>Expand Shannon Entropy documentation</strong></summary>

  A modular three-script pipeline for calculating and visualizing Shannon entropy from Multiple Sequence Alignments (MSA). Works with both amino acid (AA) and nucleotide (DNA) sequences from Read2Tree output.

## Overview

This pipeline processes consensus sequences generated by [Read2Tree](https://github.com/DessimozLab/read2tree) to calculate position-wise Shannon entropy. The analysis can be performed with or without reference sequences and supports grouping by metadata variables (e.g., genotype, subgroup, time phase).

**Pipeline workflow:**
1. `msa_to_position_table.py` - Convert MSA files to position table
2. `calculate_entropy.py` - Calculate Shannon entropy per position
3. `plot_entropy.R` - Generate publication-ready visualizations

---

## Input Files

### Required Files

#### 1. MSA Files (Read2Tree Output)
- **Location**: `MSA/[virus]/AA/` or `MSA/[virus]/DNA/`
- **Format**: Phylip format (relaxed)
- **Naming**: `OG1.fa`, `OG2.fa`, ..., `OG10.fa`
- **Source**: Consensus sequences generated by Read2Tree for each ortholog group

**Example MSA file format:**
```
 1232 634
s0252               ---SPITATV TKTRGIPSAI VCCLTGRDKY PHRGHCYILT SLTKTFMGTV
s0005               ---APITAYS QQTRGLLGCI ITSLTGRDKN QVEGEVQIVS TATQTFLATC
HepC_SRR1170677_1   ---APITAYS QQTRGLLGCI ITSLTGRDKN QVEGEVQVVS TATQSFLATC
```

**Note on sample labels:**
- Reference samples typically start with `s0` (e.g., `s0252`)
- Read samples typically start with virus prefix (e.g., `HepC_SRR1170677_1`)
- Read2Tree adds `_R1` or `_R2` suffixes (automatically cleaned by scripts)

#### 2. OG Mapping Table
- **File**: `hepc_ogs.csv` (or similar for your virus)
- **Format**: CSV with columns `OG,peptide` (or `OG,gene`)
- **Source**: Manual curation based on Read2Tree's `OG_genes.tsv`

**Example content:**
```csv
OG,peptide
OG1,NS3
OG2,NS4A
OG3,NS5B
OG4,NS5A
OG5,E2
```

**How to create:**
1. Read2Tree generates `OG_genes.tsv` during analysis
2. Manually curate to assign biological names to each OG
3. Save as CSV with `OG` and `peptide`/`gene` columns

### Optional Files

#### 3. Metadata Table (for grouping/filtering)
- **Format**: CSV with `sample_id` column
- **Purpose**: Filter samples or group entropy calculations

**Example:**
```csv
sample_id,genotype,collection_year,location
HepC_SRR1170677,GT1,2011,USA
HepC_SRR5122806,GT4,2009,Netherlands
s0252,GT7,2014,Reference
```

---

## Installation

### Requirements

**Python (≥3.7):**
```bash
pip install biopython pandas numpy
```

**R (≥4.0):**
- Run scripts once - they will prompt to install missing packages (`tidyverse`, `RColorBrewer`)
- Or install manually in R:
```r
install.packages(c("tidyverse", "RColorBrewer"))
```

---

## Usage

### Script 1: Convert MSA to Position Table

**Purpose:** Parse MSA files and create a long-format position table.

#### Basic Usage (Include All Samples)

**Amino acids:**
```bash
python msa_to_position_table.py \
    --msa_dir MSA/hepc/AA \
    --og_table hepc_ogs.csv \
    --output hepc_aa_positions.csv \
    --seq_type AA
```

**DNA:**
```bash
python msa_to_position_table.py \
    --msa_dir MSA/hepc/DNA \
    --og_table hepc_ogs.csv \
    --output hepc_dna_positions.csv \
    --seq_type DNA
```

#### Exclude Reference Sequences

To exclude reference strains (samples starting with `s0`):

```bash
python msa_to_position_table.py \
    --msa_dir MSA/hepc/AA \
    --og_table hepc_ogs.csv \
    --output hepc_aa_no_refs.csv \
    --seq_type AA \
    --exclude_pattern s0
```

#### Filter by Metadata

Process only samples from a specific genotype:

```bash
python msa_to_position_table.py \
    --msa_dir MSA/hepc/AA \
    --og_table hepc_ogs.csv \
    --output hepc_aa_gt1_only.csv \
    --seq_type AA \
    --metadata metadata.csv \
    --filter_column genotype \
    --filter_value GT1
```

#### All Options

```bash
python msa_to_position_table.py \
    --msa_dir MSA/hepc/AA          # Directory with OG*.fa files
    --og_table hepc_ogs.csv         # OG to gene mapping
    --output positions.csv          # Output file
    --seq_type AA                   # AA or DNA
    --exclude_pattern s0            # Pattern to exclude (optional)
    --include_all                   # Override exclusion (optional)
    --metadata metadata.csv         # Metadata for filtering (optional)
    --filter_column genotype        # Column to filter (optional)
    --filter_value GT1              # Value to keep (optional)
```

**Output:** CSV with columns:
- `sample_id` - Sample identifier
- `position` - Position in alignment (1-indexed)
- `character` - Amino acid or nucleotide at this position
- `og` - Ortholog group (e.g., OG1)
- `gene` - Gene/protein name (e.g., NS3)
- `seq_type` - AA or DNA

---

### Script 2: Calculate Shannon Entropy

**Purpose:** Calculate Shannon entropy for each position.

#### Basic Usage (No Grouping)

```bash
python calculate_entropy.py \
    --input hepc_aa_positions.csv \
    --output hepc_aa_entropy.csv
```

#### Group by Genotype

Calculate separate entropy values for each genotype:

```bash
python calculate_entropy.py \
    --input hepc_aa_positions.csv \
    --output hepc_aa_entropy_by_genotype.csv \
    --metadata metadata.csv \
    --group_by genotype
```

#### Group by Multiple Variables

For example, subgroup and time phase (useful for RSV):

```bash
python calculate_entropy.py \
    --input rsv_aa_positions.csv \
    --output rsv_entropy_by_subgroup_phase.csv \
    --metadata metadata.csv \
    --group_by subgroup time_phase
```

#### Advanced Options

```bash
python calculate_entropy.py \
    --input positions.csv \
    --output entropy.csv \
    --metadata metadata.csv          # Metadata for grouping
    --group_by genotype              # Group by column(s)
    --min_samples 10                 # Min samples per position (default: 5)
    --exclude_gaps                   # Exclude gaps from calculation
```

**Output:** CSV with columns:
- `gene` - Gene/protein name
- `og` - Ortholog group
- `position` - Position in alignment
- `seq_type` - AA or DNA
- `entropy` - Shannon entropy (bits)
- `n_samples` - Number of samples at this position
- `n_unique_chars` - Number of unique characters
- `gap_percent` - Percentage of gaps
- `most_common_char` - Most frequent character
- Plus any grouping columns if specified

---

### Script 3: Plot Entropy

**Purpose:** Generate publication-ready entropy plots.

#### Basic Usage

```bash
Rscript plot_entropy.R hepc_aa_entropy.csv plots/entropy
```

With explicit sequence type:
```bash
Rscript plot_entropy.R hepc_aa_entropy.csv plots/entropy AA
```

#### Output Files

The script generates:
1. **`entropy_all_genes.png`** - Combined faceted plot showing all genes
2. **`entropy_[gene].png`** - Individual plot for each gene (e.g., `entropy_NS3.png`)

**Plot features:**
- Auto-detects grouping variables (genotype, subgroup, time_phase)
- Colors lines by group if grouping detected
- Separate facets for each gene in combined plot
- High-resolution (300 dpi) PNG files

---

## Complete Workflow Examples

### Example 1: HepC Analysis (AA, No References)

```bash
# Step 1: Create position table (exclude references)
python msa_to_position_table.py \
    --msa_dir MSA/hepc/AA \
    --og_table hepc_ogs.csv \
    --output hepc_aa_positions_no_refs.csv \
    --seq_type AA \
    --exclude_pattern s0

# Step 2: Calculate entropy
python calculate_entropy.py \
    --input hepc_aa_positions_no_refs.csv \
    --output hepc_aa_entropy.csv

# Step 3: Plot
Rscript plot_entropy.R hepc_aa_entropy.csv plots/hepc_entropy AA
```

### Example 2: HepC by Genotype (DNA, With References)

```bash
# Step 1: Create position table (include all samples)
python msa_to_position_table.py \
    --msa_dir MSA/hepc/DNA \
    --og_table hepc_ogs.csv \
    --output hepc_dna_positions_all.csv \
    --seq_type DNA

# Step 2: Calculate entropy grouped by genotype
python calculate_entropy.py \
    --input hepc_dna_positions_all.csv \
    --output hepc_dna_entropy_by_genotype.csv \
    --metadata metadata.csv \
    --group_by genotype

# Step 3: Plot (will show separate lines per genotype)
Rscript plot_entropy.R hepc_dna_entropy_by_genotype.csv plots/hepc_dna_by_gt DNA
```

### Example 3: RSV Analysis (AA, By Subgroup and Phase)

```bash
# Step 1: Position table (no references)
python msa_to_position_table.py \
    --msa_dir MSA/rsv/AA \
    --og_table rsv_ogs.csv \
    --output rsv_aa_positions.csv \
    --seq_type AA \
    --exclude_pattern s0

# Step 2: Entropy by subgroup and time phase
python calculate_entropy.py \
    --input rsv_aa_positions.csv \
    --output rsv_entropy_by_subgroup_phase.csv \
    --metadata rsv_metadata.csv \
    --group_by subgroup time_phase

# Step 3: Plot
Rscript plot_entropy.R rsv_entropy_by_subgroup_phase.csv plots/rsv_entropy AA
```

---

## Understanding Shannon Entropy

**Shannon entropy** measures the uncertainty or variability at each position in an alignment:

- **0 bits**: All samples have the same character (conserved position)
- **1 bit**: Two characters equally distributed (50/50)
- **2 bits**: Four characters equally distributed (25/25/25/25)
- **Higher values**: More variability

**For amino acids:** Maximum entropy = log₂(20) ≈ 4.32 bits
**For DNA:** Maximum entropy = log₂(4) = 2 bits

---

## Customization

### Adding Domain Annotations (Advanced)

For genes with known functional domains (like RSV's G protein), you can modify `plot_entropy.R` to add domain shading:

```r
# Create domain dataframe
g_domains <- data.frame(
  aa_start = c(66, 164, 183),
  aa_stop = c(163, 182, 237),
  domain = c("Mucin-like", "Central conserved", "Hinge")
)

# Pass to plotting function
plot_entropy_per_gene(entropy_df, "G", "plots/G_protein.png", 
                      domain_df = g_domains)
```

---

## Troubleshooting

### Common Issues

**Issue:** "No package called 'tidyverse'"
```bash
# Solution: Install packages in command-line R
Rscript -e "install.packages(c('tidyverse', 'RColorBrewer'), repos='https://cran.rstudio.com/')"
```

**Issue:** "object 'entropy' not found" when plotting
- **Cause:** Trying to plot position table instead of entropy table
- **Solution:** Run `calculate_entropy.py` first (Step 2)

**Issue:** No data extracted from MSA files
- **Check:** MSA files are in correct format (Phylip relaxed)
- **Check:** OG names in files match those in `hepc_ogs.csv`

**Issue:** All samples excluded by filter
- **Check:** Pattern in `--exclude_pattern` matches your reference IDs
- **Try:** Run without `--exclude_pattern` to see all sample IDs

---

## Output Interpretation

### Position Table Example
```csv
sample_id,position,character,og,gene,seq_type
HepC_SRR1170677,1,-,OG1,NS3,AA
HepC_SRR1170677,2,-,OG1,NS3,AA
HepC_SRR1170677,3,-,OG1,NS3,AA
HepC_SRR1170677,4,A,OG1,NS3,AA
```

### Entropy Table Example
```csv
gene,og,position,seq_type,entropy,n_samples,n_unique_chars,gap_percent,most_common_char
NS3,OG1,1,AA,0.000,1200,1,100.0,-
NS3,OG1,4,AA,0.245,1200,3,0.0,A
NS3,OG1,5,AA,1.891,1200,8,0.0,P
```

---

## Citation

If you use these scripts, please cite:

- **Read2Tree**: Dylus et al. (2024) "Inference of phylogenetic trees directly from raw sequencing reads using Read2Tree" *Nature Biotechnology*
- Your study where you applied this pipeline

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
</details>
