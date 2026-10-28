# ==============================================================================
# Script 02: LULU Post-clustering Curation (Cyanobacteria 16S)
# Project: Biological Soil Crust Peru
# Purpose: Collapse redundant ASVs and generate summary indices
# ==============================================================================

# --- 1. Load Libraries ---
if (!requireNamespace("lulu", quietly = TRUE)) devtools::install_github("tobiasgf/lulu")
library(dada2)
library(lulu)
library(dplyr)
library(stringr)
library(here)
# --- 2. Define Paths ---
input_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
temp_dir   <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Cyano_LULU")

if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
setwd(temp_dir)

# --- 3. Helper Functions ---
name.change <- function(x) {
  # Fix names by stripping the fastq.gz extensions to match standard sample IDs
  names_clean <- gsub("_F_filt\\.fastq\\.gz", "", colnames(x))
  colnames(x) <- names_clean
  return(x)
}

index.info <- function(x) {
  # Calculate diversity indices (Richness and Total Reads)
  y <- data.frame(
    Sample = colnames(x),
    Richness = colSums(x > 0),
    Reads = colSums(x)
  )
  return(y)
}

# --- 4. Main LULU Processing Loop ---
strategies <- c("independent", "totalpool", "pseudopool")

for (type in strategies) {
  cat(paste0("\n======================================================\n"))
  cat(paste0("--- Running LULU for strategy: ", toupper(type), " ---\n"))
  
  # A. Load DADA2 Object
  seqtab <- readRDS(file.path(input_dir, paste0("Cyano_seqtab_nochim_", type, ".rds")))
  
  # B. Generate FASTA
  fasta_file <- file.path(temp_dir, paste0("Cyano_ASVs_", type, ".fasta"))
  uniquesToFasta(seqtab, fout = fasta_file, ids = paste0("OTU", seq(length(getSequences(seqtab)))))
  
  match_list_file <- file.path(temp_dir, paste0("Cyano_match_list_", type, ".txt"))
  
  # C. Run BLAST commands separately (Windows CMD fix)
  cat("  > Building BLAST database...\n")
  makeblastdb_cmd <- paste0("makeblastdb -in \"", fasta_file, "\" -parse_seqids -dbtype nucl")
  system(makeblastdb_cmd)
  
  cat("  > Running BLASTn alignment...\n")
  blastn_cmd <- paste0("blastn -db \"", fasta_file, "\" ",
                       "-outfmt \"6 qseqid sseqid pident\" ",
                       "-out \"", match_list_file, "\" ",
                       "-qcov_hsp_perc 80 -perc_identity 84 ",
                       "-query \"", fasta_file, "\"")
  system(blastn_cmd)
  
  # D. LULU Curation
  cat("  > Applying LULU curation algorithm...\n")
  asv_table <- as.data.frame(t(seqtab))
  rownames(asv_table) <- paste0("OTU", seq(nrow(asv_table)))
  
  match_list <- read.table(match_list_file)
  curated_result <- lulu(asv_table, match_list)
  
  # E. Reconstruct Table
  kept_indices <- as.numeric(gsub("OTU", "", rownames(curated_result$curated_table)))
  curated_table <- t(curated_result$curated_table)
  colnames(curated_table) <- colnames(seqtab)[kept_indices]
  
  # Clean names
  curated_table <- name.change(curated_table)
  
  # F. Generate and Save Indices
  indices <- index.info(curated_table)
  write.csv(indices, file.path(temp_dir, paste0("Cyano_indices_curated_", type, ".csv")), row.names = FALSE)
  
  # G. Save Curated Object
  saveRDS(curated_result, file.path(output_dir, paste0("Cyano_LULU_curated_", type, ".rds")))
  
  cat(paste0("  > Initial ASVs: ", ncol(seqtab), "\n"))
  cat(paste0("  > Curated OTUs: ", curated_result$curated_count, "\n"))
  cat(paste0("  > Discarded ASVs: ", ncol(seqtab) - curated_result$curated_count, "\n"))
}

cat("\n======================================================\n")
cat("[SUCCESS] Script 02 Completed for all strategies.\n")
cat("======================================================\n")