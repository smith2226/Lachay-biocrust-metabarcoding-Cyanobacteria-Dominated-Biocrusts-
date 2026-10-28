# ==============================================================================
# Script 05: Taxonomic Assignment (Cyanobacteria 16S - CyanoSeq V1.3)
# Project: Biological Soil Crust Peru
# Purpose: Assign taxonomy, update Master Phyloseq, and split by Ecosystem
# ==============================================================================

# --- 1. Load Libraries ---
library(dada2)
library(phyloseq)
library(Biostrings)
library(dplyr)
library(here)
# --- 2. Define Paths ---
input_dir   <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
tax_csv_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Cyano_Filtering")
tax_csv     <- file.path(tax_csv_dir, "Cyano_Taxonomy_Table.csv")
# CyanoSeq V1.3 Database Path
db_path     <- here("Zenodo_Metabarcoding", "01_Data", "Databases", "Cyanoseq", "CyanoSeqV1.3_dada2.fastq.gz")
if(!dir.exists(tax_csv_dir)) dir.create(tax_csv_dir, recursive = TRUE)
if (!file.exists(db_path)) stop("Database not found! Check the exact file name and path (ensure it is not .fasta.gz).")
# --- 3. Custom Function to Split by Ecosystem ---
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
cat("--- Running Taxonomy Assignment for Cyano (TOTALPOOL) ---\n")

# A. Load Master Phyloseq Object
ps_obj <- readRDS(file.path(input_dir, "Cyano_Phyloseq_Master_Totalpool.rds"))

# B. Extract DNA Sequences
cat("  > Extracting DNA sequences for classification...\n")
seqs <- as.character(refseq(ps_obj))

# C. Run Taxonomic Assignment (This might take a few minutes!)
cat("  > Running DADA2 assignTaxonomy using CyanoSeq V1.3...\n")
tax_info <- assignTaxonomy(seqs, db_path, multithread = TRUE, minBoot = 80)

# D. Clean Taxonomy Table
rownames(tax_info) <- taxa_names(ps_obj) # Assign "OTU1", "OTU2"... as rownames

# Standard levels for 16S databases (CyanoSeq / SILVA)
bact_ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_info) <- bact_ranks[1:ncol(tax_info)]

# E. Save Taxonomy Table as CSV
cat("  > Saving Taxonomy CSV report...\n")
write.csv(tax_info, file = tax_csv, row.names = TRUE)

# F. Merge Taxonomy into Master Phyloseq
cat("  > Integrating taxonomy into Master Phyloseq object...\n")
tax_table(ps_obj) <- tax_table(as.matrix(tax_info))

# G. Re-Split by Ecosystem
cat("  > Splitting ecosystems (now with taxonomy)...\n")
ps_ecosystem_list <- create_phyloseq_by_ecosystem(ps_obj, ecosystem_col = "Ecosystem")

# H. Save Final Objects
saveRDS(ps_obj, file.path(output_dir, "Cyano_Phyloseq_Master_Totalpool_TAX.rds"))
saveRDS(ps_ecosystem_list, file.path(output_dir, "Cyano_Phyloseq_EcosysList_Totalpool_TAX.rds"))

names(ps_ecosystem_list) <- paste0("ps_Cyano_", names(ps_ecosystem_list))
list2env(ps_ecosystem_list, envir = .GlobalEnv)

cat(paste0("\n  [SUCCESS]\n",
           "    - Taxonomy assigned for ", nrow(tax_info), " OTUs.\n",
           "    - Taxonomy table saved at: ", tax_csv, "\n",
           "    - Final TAX Objects ready for ecological analysis!\n"))
cat("======================================================\n")