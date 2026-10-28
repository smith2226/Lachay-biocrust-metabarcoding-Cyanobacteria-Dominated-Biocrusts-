# ==============================================================================
# Script 04: Phyloseq Assembly & Replicate Merging (Green Algae / G500-G800)
# Project: Biological Soil Crust Peru (FAIR)
# ==============================================================================

# --- 1. Load Libraries ---
library(phyloseq)
library(Biostrings)
library(dplyr)

# --- 1.1 Set working directory ---
setwd(normalizePath(file.path(here::here(), "..")))

# --- 2. Define Paths ---
base_dir <- normalizePath(file.path(here::here(), ".."))

input_dir  <- file.path(base_dir, "GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir <- file.path(base_dir, "GitHub_Metabarcoding", "02_Analysis_Core", "objects")
meta_path  <- file.path(base_dir, "Zenodo_Metabarcoding", "01_Data", "Metadata", "Pierr_metadata.csv")
export_dir <- file.path(base_dir, "GitHub_Metabarcoding", "02_Analysis_Core", "results", "exported_data")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

cat("\n======================================================\n")
cat("--- Building Biological-Sample Phyloseq & Exporting FASTA (Green Algae) ---\n")

# --- 3. Load Data ---
cat("  > Loading decontaminated ASV table (Totalpool)...\n")
otu_replicates <- readRDS(file.path(input_dir, "Algae_final_master_totalpool.rds"))

cat("  > Loading metadata...\n")
meta_df <- read.csv(meta_path, stringsAsFactors = FALSE)

if (!"SampleID" %in% colnames(meta_df)) {
  stop("Metadata must contain a column named 'SampleID'.")
}

rownames(meta_df) <- as.character(meta_df$SampleID)

# --- 4. Validate OTU table orientation ---
# Expected input: rows = PCR replicates, columns = ASV sequences
if (is.null(rownames(otu_replicates))) {
  stop("The OTU table has no rownames. Expected rownames like S001_rep1.")
}

if (is.null(colnames(otu_replicates))) {
  stop("The OTU table has no colnames. Expected ASV sequences as columns.")
}

otu_replicates <- as.data.frame(otu_replicates)

cat("  > Input replicate-level table dimensions:\n")
cat("    - PCR replicate rows:", nrow(otu_replicates), "\n")
cat("    - ASV sequence columns:", ncol(otu_replicates), "\n")

# --- 5. Merge PCR replicates into biological samples ---
# This matches the logic used in CYA and Euka01:
# S001_rep1 + S001_rep2 + S001_rep3 + S001_rep4 -> S001

cat("  > Merging PCR replicates into biological samples...\n")

sample_bases <- sub("_rep[1-4]$", "", rownames(otu_replicates))

if (all(sample_bases == rownames(otu_replicates))) {
  warning(
    "No '_rep1-4' suffix was detected in rownames. ",
    "The table may already be at biological-sample level."
  )
}

missing_meta <- setdiff(unique(sample_bases), rownames(meta_df))
if (length(missing_meta) > 0) {
  cat("  > SampleIDs missing in metadata:\n")
  print(missing_meta)
  stop("Some biological sample IDs inferred from replicate names are missing in metadata.")
}

merged_otu_list <- lapply(
  split(as.data.frame(otu_replicates), sample_bases),
  function(x) {
    colSums(x, na.rm = TRUE)
  }
)

otu_merged <- do.call(rbind, merged_otu_list)
otu_merged <- as.matrix(otu_merged)

# Remove ASVs with zero counts after merging
otu_merged <- otu_merged[, colSums(otu_merged) > 0, drop = FALSE]

cat("  > Biological-sample table dimensions after merging:\n")
cat("    - Biological samples:", nrow(otu_merged), "\n")
cat("    - ASVs retained:", ncol(otu_merged), "\n")
cat("    - Total reads:", sum(otu_merged), "\n")

# --- 6. Align metadata with merged OTU table ---
cat("  > Aligning metadata with biological samples...\n")

meta_filtered <- meta_df[rownames(otu_merged), , drop = FALSE]

if (!all(rownames(meta_filtered) == rownames(otu_merged))) {
  stop("Metadata rownames do not match merged OTU table rownames.")
}

# --- 7. Create Base Phyloseq Object ---
cat("  > Assembling base phyloseq object at biological-sample level...\n")

ps_base <- phyloseq(
  otu_table(otu_merged, taxa_are_rows = FALSE),
  sample_data(meta_filtered)
)

print(ps_base)

# --- 8. Rename ASVs and preserve original sequences ---
cat("  > Renaming ASVs and extracting sequences...\n")

dna_seqs <- taxa_names(ps_base)

if (is.null(dna_seqs) || length(dna_seqs) == 0) {
  stop("No ASV sequences found as taxa names.")
}

asv_names <- paste0("ASV_", seq_along(dna_seqs))

# Add reference sequences before renaming
refseq_obj <- DNAStringSet(dna_seqs)
names(refseq_obj) <- asv_names

taxa_names(ps_base) <- asv_names
ps_base <- merge_phyloseq(ps_base, refseq(refseq_obj))

# --- 9. Export FASTA for taxonomy ---
cat("  > Exporting clean FASTA file for Script 05...\n")

writeXStringSet(
  refseq_obj,
  file.path(export_dir, "Algae_ASVs_for_Taxonomy.fasta")
)

# --- 10. Save objects and sequence mapping ---
cat("  > Saving biological-sample base phyloseq and sequence mapping...\n")

saveRDS(
  ps_base,
  file.path(output_dir, "Algae_Phyloseq_Base.rds")
)

mapping_df <- data.frame(
  ASV_ID = asv_names,
  Sequence = dna_seqs,
  stringsAsFactors = FALSE
)

write.csv(
  mapping_df,
  file.path(export_dir, "Algae_ASV_Sequence_Mapping.csv"),
  row.names = FALSE
)

cat("\n[SUCCESS] Script 04 completed.\n")
cat("GreenAlgae base phyloseq is now at biological-sample level.\n")
cat("FASTA ready for PR2 taxonomy.\n")
cat("======================================================\n")
