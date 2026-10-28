# ==============================================================================
# Script 04: Phyloseq Assembly & Replicate Merging (Cyanobacteria 16S) - FIXED
# Project: Biological Soil Crust Peru
# Purpose: Merge PCR replicates, generate reference FASTA, and build Phyloseq
# ==============================================================================

# --- 1. Load Libraries ---
library(phyloseq)
library(Biostrings)
library(dplyr)
library(here)
# --- 2. Define Paths (Zenodo FAIR Structure) ---
input_dir   <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
meta_path   <- here("Zenodo_Metabarcoding", "01_Data", "Metadata", "Pierr_metadata.csv")
fasta_dir   <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "exported_data")
# Create output directories if they do not exist
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if(!dir.exists(fasta_dir)) dir.create(fasta_dir, recursive = TRUE)
# --- 3. Custom Function to Split by Ecosystem (Robust Method) ---
create_phyloseq_by_ecosystem <- function(physeq, ecosystem_col = "Ecosystem") {
  sd_df <- as(sample_data(physeq), "data.frame")
  if (ecosystem_col %in% colnames(sd_df)) {
    sd_df[[ecosystem_col]] <- trimws(sd_df[[ecosystem_col]])
    sample_data(physeq) <- sample_data(sd_df)
  } else { stop(paste("Column", ecosystem_col, "was not found in metadata.")) }
  
  ecosystems <- unique(as.character(sample_data(physeq)[[ecosystem_col]]))
  ecosystems <- ecosystems[!is.na(ecosystems) & ecosystems != ""]
  
  ps_list <- list()
  for (eco in ecosystems) {
    cat(paste("  > Creating clean subset for:", eco, "...\n"))
    samples_to_keep <- rownames(sd_df)[sd_df[[ecosystem_col]] == eco]
    ps_sub <- prune_samples(samples_to_keep, physeq)
    ps_list[[eco]] <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub) 
  }
  return(ps_list)
}

# --- 4. Main Processing ---
cat("\n======================================================\n")
cat("--- Assembling Phyloseq for Cyano (TOTALPOOL) ---\n")

# A. Load Filtered Matrix (from Script 03) & Metadata
otu_replicates <- readRDS(file.path(input_dir, "Cyano_final_master_totalpool.rds"))
exdata <- read.csv(meta_path, stringsAsFactors = FALSE)
if("SampleID" %in% colnames(exdata)){
  rownames(exdata) <- exdata$SampleID
} else { stop("Metadata must contain a column named 'SampleID'.") }

# ---- PARCHE CRÍTICO: RECUPERAR SECUENCIAS DE ADN ----
cat("  > Recovering actual DNA sequences from original DADA2 object...\n")
original_dada2 <- readRDS(file.path(input_dir, "Cyano_seqtab_nochim_totalpool.rds"))
original_seqs <- colnames(original_dada2)
names(original_seqs) <- paste0("OTU", seq_along(original_seqs))

current_otu_ids <- colnames(otu_replicates) # "OTU1", "OTU14", etc.
real_dna_sequences <- original_seqs[current_otu_ids] # Cruza los IDs con el ADN
colnames(otu_replicates) <- real_dna_sequences # Restaura el ADN a las columnas
# -----------------------------------------------------

# B. Merge PCR Replicates (colSums by Biological Sample Base Name)
cat("  > Merging PCR replicates into biological samples...\n")
sample_bases <- sub("_rep[1-4]$", "", rownames(otu_replicates))
merged_otu_list <- lapply(split(as.data.frame(otu_replicates), sample_bases), colSums)
merged_otu_matrix <- do.call(rbind, merged_otu_list)

# C. Filter Metadata to match surviving samples
surviving_samples <- rownames(merged_otu_matrix)
exdata_subset <- exdata[surviving_samples, , drop = FALSE]

# D. Generate Reference FASTA (CRITICAL for Script 05 Taxonomy)
cat("  > Extracting DNA sequences and generating FASTA file...\n")
actual_sequences <- colnames(merged_otu_matrix)
otu_clean_ids <- paste0("OTU", seq_along(actual_sequences))

refseq_obj <- Biostrings::DNAStringSet(actual_sequences)
names(refseq_obj) <- otu_clean_ids

fasta_file_path <- file.path(fasta_dir, "Cyano_refseq_totalpool.fasta")
Biostrings::writeXStringSet(refseq_obj, filepath = fasta_file_path)

# E. Build the Master Phyloseq Object
cat("  > Building the Master Phyloseq Object...\n")
colnames(merged_otu_matrix) <- otu_clean_ids

ps_obj <- phyloseq(
  otu_table(merged_otu_matrix, taxa_are_rows = FALSE), 
  sample_data(exdata_subset),
  refseq(refseq_obj)
)

# F. Split by Ecosystem 
cat("  > Splitting object by ecosystems...\n")
ps_ecosystem_list <- create_phyloseq_by_ecosystem(ps_obj, ecosystem_col = "Ecosystem")

# G. Save the Objects & Push to Environment
saveRDS(ps_obj, file.path(output_dir, "Cyano_Phyloseq_Master_Totalpool.rds"))
saveRDS(ps_ecosystem_list, file.path(output_dir, "Cyano_Phyloseq_EcosysList_Totalpool.rds"))

names(ps_ecosystem_list) <- paste0("ps_Cyano_", names(ps_ecosystem_list))
list2env(ps_ecosystem_list, envir = .GlobalEnv)

# H. Final Summary
cat(paste0("\n  [SUCCESS]\n",
           "    - Total Biological Samples: ", nsamples(ps_obj), "\n",
           "    - Total Taxa (OTUs): ", ntaxa(ps_obj), "\n",
           "    - Ecosystems generated: ", paste(names(ps_ecosystem_list), collapse = ", "), "\n",
           "    - FASTA saved at: ", fasta_file_path, "\n"))

cat("======================================================\n")