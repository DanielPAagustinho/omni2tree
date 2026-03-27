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
  -t template_v5.html \
  -l demo2 \
  -o demo2/output/omni2treeview_demo2
```

This command will create an HTML file in the specified output directory that visualizes the phylogenetic tree along with the provided metadata.

Options:
- `-n`: Path to the Newick tree file.
- `-m`: Path to the metadata CSV file.
- `-t`: Path to the HTML template file.
- `-l`: Label for the dataset.
- `-o`: Output prefix for the generated HTML file and other related files. 


## Input files:

**Newick tree file (.nwk):** output from Omni2tree.

**Metadata CSV file (.csv):** Metadata in CSV format for each sample. The first column should be `sample_id`, the second column should be `label`, and the `label` values should match the leaf names in the Newick tree file.

**Entropy folder (optional):** If a folder named `entropy` exists beside the Newick tree file, Omni2treeView will automatically copy it to the output directory and expose the PNG files in the HTML entropy panel.

**HTML template file (.html):** An HTML template file for rendering the tree view.


## Output files:

The output will include:
- An HTML file for interactive visualization (html).
- A JSON file for tree structure (tree.json).
- A JSON file for global tree settings (tree_meta.json).
- A integrated CSV file combining tree and metadata information (meta.csv).
