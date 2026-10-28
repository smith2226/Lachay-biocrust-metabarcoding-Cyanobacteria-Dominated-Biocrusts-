# ==============================================================================
# Script 06: Taxonomic Guild Subsetting (Cyanobacteria 16S)
# Project: Biological Soil Crust Peru
# Purpose: Extract True Cyanobacteria and Chloroplasts for ecological Papers
# ==============================================================================

# --- 1. Load Libraries ---
library(phyloseq)
library(dplyr)
library(here)
# --- 2. Define Paths ---
input_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects", "subsets")

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 3. Load Data ---
ps_list_file <- file.path(input_dir, "Cyano_Phyloseq_EcosysList_Totalpool_TAX.rds")
if (file.exists(ps_list_file)) {
  ps_list <- readRDS(ps_list_file)
  # Ensure they load with the exact names from Script 05 (ps_Cyano_Lomas, etc.)
  list2env(ps_list, envir = .GlobalEnv)
} else {
  stop("Master list not found. Run Script 05 first.")
}

# --- 4. Robust Functions to Extract Cyano and Chloroplasts ---
extract_true_cyano <- function(physeq) {
  if (is.null(physeq)) return(NULL)
  tax_matrix <- as(tax_table(physeq), "matrix")
  tax_df <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
  
  # Logic: Must be Phylum Cyanobacteriota AND Class MUST NOT be Chloroplast
  is_cyano <- !is.na(tax_df$Phylum) & tax_df$Phylum == "Cyanobacteriota"
  is_not_chloro <- is.na(tax_df$Class) | tax_df$Class != "Chloroplast"
  
  target_otus <- rownames(tax_df)[is_cyano & is_not_chloro]
  
  if (length(target_otus) == 0) return(NULL)
  ps_sub <- prune_taxa(target_otus, physeq)
  return(prune_taxa(taxa_sums(ps_sub) > 0, ps_sub))
}

extract_chloroplasts <- function(physeq) {
  if (is.null(physeq)) return(NULL)
  tax_matrix <- as(tax_table(physeq), "matrix")
  tax_df <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
  
  # Logic: Class MUST be Chloroplast
  is_chloro <- !is.na(tax_df$Class) & tax_df$Class == "Chloroplast"
  
  target_otus <- rownames(tax_df)[is_chloro]
  
  if (length(target_otus) == 0) return(NULL)
  ps_sub <- prune_taxa(target_otus, physeq)
  return(prune_taxa(taxa_sums(ps_sub) > 0, ps_sub))
}

# --- 5. Process Ecosystems ---
cat("\n======================================================\n")
cat("--- Splitting Cyano Guilds by Ecosystem ---\n")

# A. LOMAS
cat("\n  > Processing LOMAS...\n")
if(exists("ps_Cyano_Lomas")) {
  ps_Lomas_TrueCyano <- extract_true_cyano(ps_Cyano_Lomas)
  cat("    - True Cyanobacteria: ", ifelse(is.null(ps_Lomas_TrueCyano), 0, ntaxa(ps_Lomas_TrueCyano)), "OTUs\n")
  ps_Lomas_Chloro <- extract_chloroplasts(ps_Cyano_Lomas)
  cat("    - Chloroplasts:       ", ifelse(is.null(ps_Lomas_Chloro), 0, ntaxa(ps_Lomas_Chloro)), "OTUs\n")
  
  if(!is.null(ps_Lomas_TrueCyano)) saveRDS(ps_Lomas_TrueCyano, file.path(output_dir, "Cyano_Lomas_TrueCyano.rds"))
  if(!is.null(ps_Lomas_Chloro)) saveRDS(ps_Lomas_Chloro, file.path(output_dir, "Cyano_Lomas_Chloroplasts.rds"))
} else { cat("    [!] ps_Cyano_Lomas not found.\n") }

# B. GLACIER
cat("\n  > Processing GLACIER...\n")
if(exists("ps_Cyano_Glacier")) {
  ps_Glacier_TrueCyano <- extract_true_cyano(ps_Cyano_Glacier)
  cat("    - True Cyanobacteria: ", ifelse(is.null(ps_Glacier_TrueCyano), 0, ntaxa(ps_Glacier_TrueCyano)), "OTUs\n")
  ps_Glacier_Chloro <- extract_chloroplasts(ps_Cyano_Glacier)
  cat("    - Chloroplasts:       ", ifelse(is.null(ps_Glacier_Chloro), 0, ntaxa(ps_Glacier_Chloro)), "OTUs\n")
  
  if(!is.null(ps_Glacier_TrueCyano)) saveRDS(ps_Glacier_TrueCyano, file.path(output_dir, "Cyano_Glacier_TrueCyano.rds"))
  if(!is.null(ps_Glacier_Chloro)) saveRDS(ps_Glacier_Chloro, file.path(output_dir, "Cyano_Glacier_Chloroplasts.rds"))
} else { cat("    [!] ps_Cyano_Glacier not found.\n") }

# C. GLACIER EXPERIMENT
cat("\n  > Processing GLACIER EXPERIMENT...\n")
if(exists("ps_Cyano_Glacier_Experiment")) {
  ps_GlacierExp_TrueCyano <- extract_true_cyano(ps_Cyano_Glacier_Experiment)
  cat("    - True Cyanobacteria: ", ifelse(is.null(ps_GlacierExp_TrueCyano), 0, ntaxa(ps_GlacierExp_TrueCyano)), "OTUs\n")
  ps_GlacierExp_Chloro <- extract_chloroplasts(ps_Cyano_Glacier_Experiment)
  cat("    - Chloroplasts:       ", ifelse(is.null(ps_GlacierExp_Chloro), 0, ntaxa(ps_GlacierExp_Chloro)), "OTUs\n")
  
  if(!is.null(ps_GlacierExp_TrueCyano)) saveRDS(ps_GlacierExp_TrueCyano, file.path(output_dir, "Cyano_GlacierExp_TrueCyano.rds"))
  if(!is.null(ps_GlacierExp_Chloro)) saveRDS(ps_GlacierExp_Chloro, file.path(output_dir, "Cyano_GlacierExp_Chloroplasts.rds"))
} else { cat("    [!] ps_Cyano_Glacier_Experiment not found.\n") }

cat("\n======================================================\n")
cat("[SUCCESS] Script 06 completed! Guild subsets saved to:\n")
cat(output_dir, "\n")
cat("======================================================\n")