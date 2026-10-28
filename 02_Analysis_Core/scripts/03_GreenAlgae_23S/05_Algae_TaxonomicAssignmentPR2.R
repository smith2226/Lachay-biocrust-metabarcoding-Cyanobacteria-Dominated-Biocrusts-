# ==============================================================================
# Script 05: Taxonomic Assignment PR2 - Green Algae / G500-G800
# ==============================================library(phyloseq)
library(dada2)
library(Biostrings)
library(dplyr)

# --- Set base directory ---
setwd(normalizePath(file.path(here::here(), "..")))

base_dir <- normalizePath(file.path(here::here(), ".."))

export_dir <- file.path(
  base_dir,
  "GitHub_Metabarcoding",
  "02_Analysis_Core",
  "results",
  "exported_data"
)

output_dir <- file.path(
  base_dir,
  "GitHub_Metabarcoding",
  "02_Analysis_Core",
  "objects"
)

db_path <- file.path(
  base_dir,
  "Zenodo_Metabarcoding",
  "01_Data",
  "Databases",
  "PR2",
  "pr2_version_5.1.0_SSU_dada2.fasta.gz"
)

# --- Check files ---
cat("Checking input files...\n")
cat("FASTA:", file.path(export_dir, "Algae_ASVs_for_Taxonomy.fasta"), "\n")
cat("PR2 database:", db_path, "\n")
cat("Phyloseq base:", file.path(output_dir, "Algae_Phyloseq_Base.rds"), "\n")

if (!file.exists(file.path(export_dir, "Algae_ASVs_for_Taxonomy.fasta"))) {
  stop("ASV FASTA file not found. Run Script 04 first.")
}

if (!file.exists(file.path(output_dir, "Algae_Phyloseq_Base.rds"))) {
  stop("Algae_Phyloseq_Base.rds not found. Run Script 04 first.")
}

if (!file.exists(db_path)) {
  cat("\nPR2 database was not found at:\n")
  cat(db_path, "\n\n")
  cat("Available files in PR2 folder:\n")
  print(list.files(dirname(db_path), full.names = FALSE))
  stop("PR2 database not found. Check the exact filename.")
}
# --- 3. Load ASV Sequences ---
cat("  > Loading ASV Sequence Mapping...\n")
mapping_file <- file.path(export_dir, "Algae_ASV_Sequence_Mapping.csv")
if(!file.exists(mapping_file)) stop("Mapping file not found! Run Script 04 first.")

mapping_df <- read.csv(mapping_file, stringsAsFactors = FALSE)
# Extract real DNA sequences and name them with ASV_IDs
seqs <- mapping_df$Sequence
names(seqs) <- mapping_df$ASV_ID

cat(paste("    - Total ASVs to classify:", length(seqs), "\n"))

# --- 4. Assign Taxonomy ---
cat("  > Running assignTaxonomy against PR2 database (High CPU/RAM)...\n")
# PR2 typical levels: Domain, Supergroup, Division, Class, Order, Family, Genus, Species
taxLevels <- c("Domain", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
gc()
taxa <- assignTaxonomy(seqs, 
                       db_path, 
                       multithread = T, 
                       taxLevels = taxLevels,
                       minBoot = 50) 

# Ensure rownames are the clean ASV_IDs, not the DNA sequences
rownames(taxa) <- names(seqs)

# --- 5. Save Object ---
cat("  > Saving Taxonomy Object...\n")
saveRDS(taxa, file.path(output_dir, "Algae_tax_final.rds"))

cat("\n[SUCCESS] Script 05 Completed! Taxonomy assigned and saved.\n")
cat("======================================================\n")
