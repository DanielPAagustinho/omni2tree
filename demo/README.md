# Omni2Tree demo dataset

This demo keeps the input small while preserving enough tree structure for a useful visualization.

## Selection

- 4 readsets: `SRR24833886`, `SRR24833879`, `SRR24833888`, `SRR24833895`.
- 5 reference assemblies: `A2`, `HRSV/A/England/397/2017`, `hRSV/A/USA/202195752/2021`, `hRSV/B/Australia/VIC-RCH056/2019`, `hMPV_CAN97-83`.
- `hMPV_CAN97-83` is the outgroup.

The readsets total about 53M on disk. The selected reference AA FASTA files are all about 3-6 KB each, and the selected CDS FASTA files are about 14-16 KB each.

## Metadata columns

`data/metadata.csv` keeps columns that should be useful for visualization:

- `source`: Reference or Readset.
- `subgroup`: RSV-A, RSV-B, hMPV, or inferred readset subgroup.
- `lineage`: Goya 2024 lineage where available.
- `tree_group`: compact display group based on the full run.
- `location`, `collection_date`, `decade`.
- `bases`, `avg_spot_len`.
- `selection_note`: why the item is in the demo.

## Suggested run

Run from `data/`:

```bash

# First, let us download the readsets
o2t-sra -i reads.csv -o reads --layout SINGLE

# Second, let's create the reference database consisting on 5 assemblies
o2t-step1 -i accessions.csv -g outgroup.csv -T 3 --o2t_out ../results |& tee hRSV_step1.log

# Third, we generate the consensus sequence for each of the 4 read samples
parallel -j 4 o2t-step2 -r {1} --o2t_out ../results -T 2 ::: \
  $(ls reads/hRSV_*fastq | sort) |& tee -a hRSV_step2.log

# Fourth, let us combine the recovered sequences with the OG alignment and generate the phylogenetic tree, entropy analysis and visualization

o2t-step3 -o ../results -m metadata.csv -l hRSV_demo --seq_type aa -T 3 -r --exclude_pattern "s0" --min_samples 4 |& tee hRSV_step3.log
```

Note that, using the `--exluce_pattern` parameter, we exclude all the reference sequences from the entropy analysis.

## How to confirm success

To verify that the results are as expected, please compare your generated `results` folder—which contains the accumulated outputs from each step of the process—with the [demo/results](demo/results) folder. Pay special attention to the `OG_genes.tsv`, `aa_positions.csv`, `aa_entropy.csv`, the HTML visualization file, and the `.nwk` files in [O2T_RESULTS](demo/results/O2T_RESULTS).

### Example HTML Visualization

![O2T visualization demo](docs/demo_o2t_view.png)

