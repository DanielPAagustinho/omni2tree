# Omni2treeView

Omni2treeView is a sub-tool from Omni2tree software to visualize phylogenetic trees with associated metadata in an interactive HTML format. It supports various customization options and is designed to handle large datasets efficiently.

## Features
- Interactive visualization of phylogenetic trees.
- Integration of metadata for enhanced analysis.
- Direct use of Newick tree files and CSV metadata files from Omni2tree.

## Installation

To install Omni2treeView, clone the repository and navigate to the `view` directory:

```bash
git clone git@github.com:DanielPAagustinho/omni2tree.git
cd omni2tree/view
```

## Usage

To generate an interactive tree view, use the following command:

```bash
python3 ./omni2treeview.py \
  -n demo2/data/concat_merge_view_aa.phy.treefile \
  -m demo2/data/metadata_o2t_view.csv \
  -e demo2/data/entropy \
  -t template_v5.html \
  -l demo2 \
  -o demo2/output/omni2treeview_demo2
```

This command will create an HTML file in the specified output directory that visualizes the phylogenetic tree along with the provided metadata.

Options:
- `-n`: Path to the Newick tree file.
- `-m`: Path to the metadata CSV file.
- `-e`: Path to an entropy image directory. Optional.
- `-t`: Path to the HTML template file.
- `-l`: Tree name shown in the HTML output.
- `-o`: Output prefix for the generated files.


## Input files:

**Newick tree file (.nwk / .treefile):** Newick tree output from Omni2tree.

**Metadata CSV file (.csv):** Metadata in CSV format for each sample.

Expected format:
- Row 1: column names.
- Row 2: column types. Allowed values are `character`, `date`, `numeric`, and `integer`.
- Row 3 onward: metadata values.

Required rules:
- Column 1 must be `sample_id`.
- Column 2 must be `label`.
- `sample_id` values must be unique.
- `label` values must be unique.
- Empty values are converted to `NA`.

Matching behavior:
- Omni2treeView first tries to match each tree leaf name to metadata column 1: `sample_id`.
- If that does not match, it then tries metadata column 2: `label`.
- If neither matches, the sample metadata is filled with `NA`.

Example:

```csv
sample_id,label,Host,Host_Group,Host_Type,Collection_State,Year
character,character,character,character,character,character,date
432,432,Bovine,Bovine,Mammalian,Ohio,2024
3101271_749133,3101271_749133,Snow goose,Goose,Avian,Kansas,2023
```

**Entropy folder (optional):**
- If `-e` is provided, that directory is used.
- If `-e` is not provided, Omni2treeView looks for a folder named `entropy` beside the Newick tree file.
- PNG files in that folder are copied to the output directory and shown in the HTML entropy panel.

**HTML template file (.html):** An HTML template file for rendering the tree view.


## Output files:

The output will include:
- An interactive HTML file: `*.html`
- A tree JSON file: `*.tree.json`
- A tree metadata JSON file: `*.tree_meta.json`
- An integrated CSV file combining tree structure and metadata: `*.meta.csv`
- A ZIP package containing the generated outputs: `*.zip`
- If entropy images are available, a copied entropy folder: `*_entropy/`

Notes on the generated view:
- The clicked-node metadata panel includes `sample_id`.
- The sample details table shows `sample_id` first, followed by `Label`, then the remaining metadata columns.
