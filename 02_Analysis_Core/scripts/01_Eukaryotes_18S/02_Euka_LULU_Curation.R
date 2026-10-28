# ==============================================================================
# Script 02: LULU Post-clustering Curation (Eukaryotes 18S)
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
# --- 2. Define Paths (Corrected for Zenodo Export Structure) ---

# Where the files from Script 01 are currently located
input_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
# Where you want the LULU curated files to go (Keeping consistency with your structure)
output_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
# Where the quality reports and temporary BLAST files will be stored
temp_dir   <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Euka_LULU")

if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
setwd(temp_dir)

# --- 3. Custom Functions ---

# Function to clean sample names (based on your S001_rep1 format)
name.change <- function(x){
  # Adapt this if your names have extra suffixes from the filtering step
  rownames(x) <- gsub("_F_filt.fastq.gz", "", rownames(x))
  return(x)
}

# Function to generate summary statistics
index.info <- function(x){
  y <- data.frame(matrix(NA, nrow = nrow(x), ncol = 0))
  y$rep <- sub(".*_rep", "", rownames(x)) 
  y$sample <- sub("_rep[1-4]$", "", rownames(x))
  y$full <- rownames(x)
  y$totseq <- rowSums(x)
  y$OTUs <- rowSums(x > 0)
  return(y)
}

# --- 4. Main Processing Loop (NP, PP, psPP) ---
# Note: Ensure these suffixes match exactly what you saved in Script 01
pooling_types <- c("independent", "totalpool", "pseudopool")

for (type in pooling_types) {
  
  cat("\n--- Processing Strategy:", type, "---\n")
  
  # A. Load Data
  input_file <- file.path(input_dir, paste0("Euka_seqtab_nochim_", type, ".rds"))
  if (!file.exists(input_file)) {
    warning(paste("File not found:", input_file))
    next
  }
  seqtab <- readRDS(input_file)
  
  # B. BLAST Internal Matching (Windows Compatibility Fix)
  fasta_file <- file.path(temp_dir, paste0("Euka_", type, "_refs.fasta"))
  uniquesToFasta(seqtab, fout = fasta_file, ids = paste0("OTU", seq(ncol(seqtab))))
  
  match_list_file <- file.path(temp_dir, paste0(type, "_match_list.txt"))
  
  # Running BLAST 
  # 1. Create DB
  system(paste0("makeblastdb -in \"", fasta_file, "\" -parse_seqids -dbtype nucl"))
  
  # 2. Run Blast with escaped quotes for Windows fmt
  blast_cmd <- paste0("blastn -db \"", fasta_file, "\" ",
                      "-outfmt \"6 qseqid sseqid pident\" ",
                      "-out \"", match_list_file, "\" ",
                      "-qcov_hsp_perc 80 -perc_identity 84 ",
                      "-query \"", fasta_file, "\"")
  
  system(blast_cmd)
  # C. LULU Curation
  # Prepare table: OTUs as rows, Samples as columns
  asv_table <- as.data.frame(t(seqtab))
  rownames(asv_table) <- paste0("OTU", seq(nrow(asv_table)))
  
  match_list <- read.table(match_list_file)
  curated_result <- lulu(asv_table, match_list)
  
  # D. Reconstruct Table
  # Extract which OTUs were kept
  kept_indices <- as.numeric(gsub("OTU", "", rownames(curated_result$curated_table)))
  curated_table <- t(curated_result$curated_table)
  colnames(curated_table) <- colnames(seqtab)[kept_indices]
  
  # Clean names
  curated_table <- name.change(curated_table)
  
  # E. Generate and Save Indices
  indices <- index.info(curated_table)
  write.csv(indices, file.path(temp_dir, paste0("Euka_indices_curated_", type, ".csv")), row.names = FALSE)
  
  # F. Save Curated Object
  saveRDS(curated_table, file.path(output_dir, paste0("Euka_seqtab_curated_", type, ".rds")))
  
  # Statistics Printout
  cat(paste0("Curation Complete for ", type, ":\n"))
  cat(paste0("  OTUs before: ", ncol(seqtab), "\n"))
  cat(paste0("  OTUs after:  ", curated_result$curated_count, "\n"))
}

# --- 5. Final Workspace Save ---
save.image(file.path(output_dir, "Euka_Script02_LULU_final.RData"))
cat("\n--- Script 02: LULU Curation Complete ---")
