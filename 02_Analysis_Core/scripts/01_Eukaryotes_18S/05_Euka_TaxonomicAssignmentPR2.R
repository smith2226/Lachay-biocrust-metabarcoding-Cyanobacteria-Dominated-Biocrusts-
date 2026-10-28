# ==============================================================================
# Script 05: Taxonomic Assignment with PR2 (18S Eukaryotes)
# Project: Biological Soil Crust Peru
# Purpose: Assign taxonomy, update Master Phyloseq, and split for Papers
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
tax_csv_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Euka_Filtering")
tax_csv     <- file.path(tax_csv_dir, "Euka_Taxonomy_Table.csv")
# PR2 Database Path (Version 5.1.0)
pr2_db_path <- here("Zenodo_Metabarcoding", "01_Data", "Databases", "PR2", "pr2_version_5.1.0_SSU_dada2.fasta.gz")

if(!dir.exists(tax_csv_dir)) dir.create(tax_csv_dir, recursive = TRUE)
if (!file.exists(pr2_db_path)) stop("PR2 database not found! Check the exact file name and path in Zenodo_Metabarcoding.")

# --- 3. Custom Function to Split by Ecosystem (Robust Method) ---
create_phyloseq_by_ecosystem <- function(physeq, ecosystem_col = "Ecosystem") {
  
  # A. Clean hidden whitespaces in the metadata
  sd_df <- as(sample_data(physeq), "data.frame")
  if (ecosystem_col %in% colnames(sd_df)) {
    sd_df[[ecosystem_col]] <- trimws(sd_df[[ecosystem_col]])
    sample_data(physeq) <- sample_data(sd_df)
  } else {
    stop(paste("Column", ecosystem_col, "was not found in metadata."))
  }
  
  # B. Identify unique ecosystems
  ecosystems <- unique(as.character(sample_data(physeq)[[ecosystem_col]]))
  ecosystems <- ecosystems[!is.na(ecosystems) & ecosystems != ""]
  
  # C. Safe splitting loop
  ps_list <- list()
  for (eco in ecosystems) {
    samples_to_keep <- rownames(sd_df)[sd_df[[ecosystem_col]] == eco]
    ps_sub <- prune_samples(samples_to_keep, physeq)
    ps_list[[eco]] <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub) 
  }
  return(ps_list)
}

# --- 4. Main Processing ---
cat("\n======================================================\n")
cat("--- Starting Taxonomic Assignment (PR2 Database) ---\n")

# A. Load Master Phyloseq Object (From Script 04)
ps_obj <- readRDS(file.path(input_dir, "Euka_Phyloseq_Master_Totalpool.rds"))

# B. Extract DNA Sequences and OTU names
cat("  > Extracting sequences for DADA2...\n")
seqs <- as.character(refseq(ps_obj))
otu_names <- names(seqs) # "OTU1", "OTU2", etc.

# C. Run Taxonomic Assignment
cat("  > Running taxonomic assignment with PR2... (This step may take 15 to 60 minutes)\n")
# Using minBoot = 80 (Strict Q1 standard for reliable assignments)
tax_info <- assignTaxonomy(seqs, pr2_db_path, multithread = TRUE, minBoot = 80)

# D. Clean Taxonomy Table and format for Phyloseq
cat("  > Formateando tabla de taxonomía...\n")
rownames(tax_info) <- otu_names # Asignamos "OTU1", "OTU2" como rownames

# Lista maestra de niveles de PR2 (hasta 9 niveles)
pr2_ranks <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

# Asignamos exactamente la cantidad de nombres que tu tabla necesita
colnames(tax_info) <- pr2_ranks[1:ncol(tax_info)]

# E. Save Taxonomy Table as CSV (Material Suplementario para tu paper)
write.csv(tax_info, file = tax_csv, row.names = TRUE)

# F. Merge Taxonomy into Master Phyloseq
cat("  > Integrando taxonomía al objeto Phyloseq...\n")
tax_table(ps_obj) <- tax_table(as.matrix(tax_info))

# G. Re-Split by Ecosystem (Now with Taxonomy included)
cat("  > Dividiendo los ecosistemas (con taxonomía incluida)...\n")
ps_ecosystem_list <- create_phyloseq_by_ecosystem(ps_obj, ecosystem_col = "Ecosystem")

# H. Save Final Objects
saveRDS(ps_obj, file.path(output_dir, "Euka_Phyloseq_Master_Totalpool_TAX.rds"))
saveRDS(ps_ecosystem_list, file.path(output_dir, "Euka_Phyloseq_EcosysList_Totalpool_TAX.rds"))

# Push to Global Environment
names(ps_ecosystem_list) <- paste0("ps_", names(ps_ecosystem_list))
list2env(ps_ecosystem_list, envir = .GlobalEnv)

cat(paste0("\n  [ASIGNACIÓN EXITOSA]\n",
           "    - Taxonomía asignada a ", nrow(tax_info), " OTUs.\n",
           "    - Tabla guardada en: ", tax_csv, "\n",
           "    - Ecosistemas actualizados listos en memoria.\n"))
cat("======================================================\n")
