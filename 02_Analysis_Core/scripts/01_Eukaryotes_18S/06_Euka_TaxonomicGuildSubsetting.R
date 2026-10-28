# ==============================================================================
# Script 06: Taxonomic Guild Subsetting (Eukaryotes 18S)
# Project: Biological Soil Crust Peru
# Purpose: Extract Fungi, Algae, and Metazoa for specific ecological Papers
# ==============================================================================

# --- 1. Load Libraries ---
library(phyloseq)
library(dplyr)
library(here)
# --- 2. Define Paths ---
input_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects", "subsets")

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 3. Load Data (If not already in Environment) ---
# If ps_Lomas, ps_Glacier, and ps_Glacier_Experiment are missing from your Environment,
# uncomment the next two lines to load them from your hard drive:
# ps_list <- readRDS(file.path(input_dir, "Euka_Phyloseq_EcosysList_Totalpool_TAX.rds"))
# list2env(ps_list, envir = .GlobalEnv)

# --- 4. Robust Function to Extract Guilds ---
extract_guild <- function(physeq, rank_name, target_value) {
  # Extract the taxonomy table as a matrix
  tax_matrix <- as(tax_table(physeq), "matrix")
  
  # Identify OTUs matching the target (using which() safely handles NAs in taxonomy)
  target_otus <- rownames(tax_matrix)[which(tax_matrix[, rank_name] == target_value)]
  
  if (length(target_otus) == 0) {
    cat(paste("    [!] No taxa found for", target_value, "in", rank_name, "\n"))
    return(NULL)
  }
  
  # Prune the phyloseq object to keep only the target OTUs
  ps_sub <- prune_taxa(target_otus, physeq)
  
  # Remove OTUs that have 0 total abundance in this specific subset
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  return(ps_sub)
}

# --- 5. Main Subsetting Process ---
cat("\n======================================================\n")
cat("--- Subsetting Taxonomic Guilds for Publications ---\n")

# ---------------------------------------------------------
# A. LOMAS ECOSYSTEM (Papers 01 & 02)
# ---------------------------------------------------------
cat("\n  > Processing LOMAS Ecosystem...\n")

ps_Lomas_Fungi <- extract_guild(ps_Lomas, "Subdivision", "Fungi")
cat("    - Fungi:      ", ifelse(is.null(ps_Lomas_Fungi), 0, ntaxa(ps_Lomas_Fungi)), "OTUs\n")

ps_Lomas_Algae <- extract_guild(ps_Lomas, "Division", "Chlorophyta")
cat("    - Microalgae: ", ifelse(is.null(ps_Lomas_Algae), 0, ntaxa(ps_Lomas_Algae)), "OTUs\n")

ps_Lomas_Fauna <- extract_guild(ps_Lomas, "Subdivision", "Metazoa")
cat("    - Microfauna: ", ifelse(is.null(ps_Lomas_Fauna), 0, ntaxa(ps_Lomas_Fauna)), "OTUs\n")

if(!is.null(ps_Lomas_Fungi)) saveRDS(ps_Lomas_Fungi, file.path(output_dir, "Euka_Lomas_Fungi.rds"))
if(!is.null(ps_Lomas_Algae)) saveRDS(ps_Lomas_Algae, file.path(output_dir, "Euka_Lomas_Algae.rds"))
if(!is.null(ps_Lomas_Fauna)) saveRDS(ps_Lomas_Fauna, file.path(output_dir, "Euka_Lomas_Fauna.rds"))

# ---------------------------------------------------------
# B. GLACIER CHRONOSEQUENCE (Paper 03)
# ---------------------------------------------------------
cat("\n  > Processing GLACIER Ecosystem...\n")

ps_Glacier_Fungi <- extract_guild(ps_Glacier, "Subdivision", "Fungi")
cat("    - Fungi:      ", ifelse(is.null(ps_Glacier_Fungi), 0, ntaxa(ps_Glacier_Fungi)), "OTUs\n")

ps_Glacier_Algae <- extract_guild(ps_Glacier, "Division", "Chlorophyta")
cat("    - Microalgae: ", ifelse(is.null(ps_Glacier_Algae), 0, ntaxa(ps_Glacier_Algae)), "OTUs\n")

ps_Glacier_Fauna <- extract_guild(ps_Glacier, "Subdivision", "Metazoa")
cat("    - Microfauna: ", ifelse(is.null(ps_Glacier_Fauna), 0, ntaxa(ps_Glacier_Fauna)), "OTUs\n")

if(!is.null(ps_Glacier_Fungi)) saveRDS(ps_Glacier_Fungi, file.path(output_dir, "Euka_Glacier_Fungi.rds"))
if(!is.null(ps_Glacier_Algae)) saveRDS(ps_Glacier_Algae, file.path(output_dir, "Euka_Glacier_Algae.rds"))
if(!is.null(ps_Glacier_Fauna)) saveRDS(ps_Glacier_Fauna, file.path(output_dir, "Euka_Glacier_Fauna.rds"))

# ---------------------------------------------------------
# C. GLACIER EXPERIMENT (Paper 04)
# ---------------------------------------------------------
cat("\n  > Processing GLACIER EXPERIMENT...\n")

ps_GlacierExp_Fungi <- extract_guild(ps_Glacier_Experiment, "Subdivision", "Fungi")
cat("    - Fungi:      ", ifelse(is.null(ps_GlacierExp_Fungi), 0, ntaxa(ps_GlacierExp_Fungi)), "OTUs\n")

ps_GlacierExp_Algae <- extract_guild(ps_Glacier_Experiment, "Division", "Chlorophyta")
cat("    - Microalgae: ", ifelse(is.null(ps_GlacierExp_Algae), 0, ntaxa(ps_GlacierExp_Algae)), "OTUs\n")

ps_GlacierExp_Fauna <- extract_guild(ps_Glacier_Experiment, "Subdivision", "Metazoa")
cat("    - Microfauna: ", ifelse(is.null(ps_GlacierExp_Fauna), 0, ntaxa(ps_GlacierExp_Fauna)), "OTUs\n")

if(!is.null(ps_GlacierExp_Fungi)) saveRDS(ps_GlacierExp_Fungi, file.path(output_dir, "Euka_GlacierExp_Fungi.rds"))
if(!is.null(ps_GlacierExp_Algae)) saveRDS(ps_GlacierExp_Algae, file.path(output_dir, "Euka_GlacierExp_Algae.rds"))
if(!is.null(ps_GlacierExp_Fauna)) saveRDS(ps_GlacierExp_Fauna, file.path(output_dir, "Euka_GlacierExp_Fauna.rds"))

cat("\n[SUCCESS] Taxonomic subsets successfully generated and saved in the 'Subsets' folder!\n")
cat("======================================================\n")