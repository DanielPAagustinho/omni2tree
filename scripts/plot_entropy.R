#!/usr/bin/env Rscript
# Plot Shannon entropy from entropy table
# Usage: Rscript plot_entropy.R entropy_table.csv output_directory

# Function to check and install required packages
check_and_install_packages <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  if (length(missing_packages) > 0) {
    cat("\n=== Missing Required Packages ===\n")
    cat("The following packages are required but not installed:\n")
    cat(paste("  -", missing_packages, collapse = "\n"), "\n\n")
    
    # Check if running interactively or from command line
    if (interactive()) {
      response <- readline(prompt = "Install these packages now? (y/n): ")
    } else {
      cat("Install these packages now? (y/n): ")
      response <- readLines("stdin", n = 1)
    }
    
    if (tolower(response) == "y") {
      cat("\nInstalling packages...\n")
      install.packages(missing_packages, repos = "https://cran.rstudio.com/", 
                       quiet = FALSE, dependencies = TRUE)
      cat("Installation complete!\n\n")
    } else {
      cat("\nPackage installation cancelled.\n")
      cat("Please install packages manually with:\n")
      cat(paste0("  install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n"))
      quit(status = 1)
    }
  }
}

# Check for required packages
required_packages <- c("tidyverse", "RColorBrewer")
check_and_install_packages(required_packages)

# Load libraries
library(tidyverse)
library(RColorBrewer)

fixed_entropy_cols <- c(
  "gene", "og", "position", "seq_type", "entropy",
  "n_samples", "n_unique_chars", "gap_percent", "most_common_char"
)

detect_group_cols <- function(df) {
  setdiff(colnames(df), fixed_entropy_cols)
}

# Function to create entropy plot for all genes
plot_entropy_all_genes <- function(entropy_df, output_file, title = "Shannon Entropy Across Positions") {
  
  group_cols <- detect_group_cols(entropy_df)
  has_groups <- length(group_cols) > 0
  
  if (has_groups) {
    group_col <- group_cols[1]
    
    p <- ggplot(entropy_df, aes(x = position, y = entropy, color = .data[[group_col]])) +
      geom_line(aes(group = .data[[group_col]])) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_brewer(palette = "Set1", name = str_to_title(gsub("_", " ", group_col))) +
      facet_wrap(~gene, scales = "free_x", ncol = 3)
  } else {
    p <- ggplot(entropy_df, aes(x = position, y = entropy)) +
      geom_line(color = "steelblue") +
      geom_point(size = 0.5, alpha = 0.7) +
      facet_wrap(~gene, scales = "free_x", ncol = 3)
  }
  
  p <- p +
    labs(
      title = title,
      x = "Position",
      y = "Entropy (bits)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  ggsave(output_file, p, width = 14, height = 10, dpi = 300, bg = "white")
  
  return(p)
}

# Function to create individual gene plots
plot_entropy_per_gene <- function(entropy_df, gene_name, output_file, 
                                   domain_df = NULL) {
  
  gene_data <- entropy_df %>% filter(gene == gene_name)
  
  if (nrow(gene_data) == 0) {
    warning(paste("No data found for gene:", gene_name))
    return(NULL)
  }
  
  group_cols <- detect_group_cols(gene_data)
  has_groups <- length(group_cols) > 0
  
  p <- ggplot(gene_data, aes(x = position, y = entropy))
  
  # Add domain annotations if provided
  if (!is.null(domain_df)) {
    p <- p + 
      geom_rect(
        data = domain_df,
        aes(xmin = aa_start, xmax = aa_stop, fill = domain),
        ymin = 0, ymax = Inf, alpha = 0.2, inherit.aes = FALSE
      ) +
      scale_fill_brewer(palette = "Accent", name = "Domain")
  }
  
  if (has_groups) {
    group_col <- group_cols[1]
    p <- p +
      geom_line(aes(color = .data[[group_col]], group = .data[[group_col]])) +
      geom_point(aes(color = .data[[group_col]]), size = 1) +
      scale_color_brewer(palette = "Set1", name = str_to_title(gsub("_", " ", group_col)))
  } else {
    p <- p +
      geom_line(color = "steelblue") +
      geom_point(size = 1, alpha = 0.7)
  }
  
  p <- p +
    labs(
      title = paste("Shannon Entropy -", gene_name, "protein"),
      x = paste(gene_name, "position"),
      y = "Entropy (bits)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  ggsave(output_file, p, width = 10, height = 6, dpi = 300, bg = "white")
  
  return(p)
}

# Main script
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript plot_entropy.R <entropy_table.csv> <output_directory> [seq_type]\n")
    cat("\nArguments:\n")
    cat("  entropy_table.csv  : CSV file from calculate_entropy.py\n")
    cat("  output_directory   : Directory to save plots\n")
    cat("  seq_type          : Optional, 'AA' or 'DNA' (default: auto-detect)\n")
    cat("\nExample:\n")
    cat("  Rscript plot_entropy.R hepc_entropy.csv plots/entropy AA\n")
    quit(status = 1)
  }
  
  entropy_file <- args[1]
  output_dir <- args[2]
  seq_type <- if (length(args) >= 3) args[3] else NULL
  
  # Check input file exists
  if (!file.exists(entropy_file)) {
    stop(paste("Error: File not found:", entropy_file))
  }
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load entropy data
  cat("Loading entropy data from:", entropy_file, "\n")
  entropy_df <- read_csv(entropy_file, show_col_types = FALSE)
  
  # Auto-detect sequence type if not provided
  if (is.null(seq_type) && "seq_type" %in% colnames(entropy_df)) {
    seq_type <- unique(entropy_df$seq_type)[1]
  }
  
  cat("Sequence type:", seq_type, "\n")
  cat("Genes found:", paste(unique(entropy_df$gene), collapse = ", "), "\n")
  
  # Check for grouping
  group_cols <- detect_group_cols(entropy_df)
  if (length(group_cols) > 0) {
    cat("Grouping by:", paste(group_cols, collapse = ", "), "\n")
  }
  
  # Plot all genes together
  cat("\nGenerating combined plot for all genes...\n")
  all_genes_plot <- file.path(output_dir, "entropy_all_genes.png")
  plot_entropy_all_genes(entropy_df, all_genes_plot)
  cat("  Saved:", all_genes_plot, "\n")
  
  # Plot individual genes
  cat("\nGenerating individual gene plots...\n")
  genes <- unique(entropy_df$gene)
  
  for (gene in genes) {
    gene_file <- file.path(output_dir, paste0("entropy_", gene, ".png"))
    plot_entropy_per_gene(entropy_df, gene, gene_file)
    cat("  Saved:", gene_file, "\n")
  }
  
  cat("\nAll plots saved to:", output_dir, "\n")
  cat("Done!\n")
}

# Run main function
if (!interactive()) {
  main()
}
