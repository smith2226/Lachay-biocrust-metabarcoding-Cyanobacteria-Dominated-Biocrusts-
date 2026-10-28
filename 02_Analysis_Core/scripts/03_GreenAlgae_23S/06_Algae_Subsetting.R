# ==============================================================================
# Script 06: Taxonomic Integration & Ecosystem Subsetting (Green Algae / G500-G800)
# Project: Biological Soil Crust Peru (FAIR)
# Purpose:
#   Integrate PR2 taxonomy into the biological-sample phyloseq object,
#   retain target green-algal plastid lineages, and generate ecosystem subsets.
# ==============================================================================

# --- 1. Load Libraries ---
library(phyloseq)
library(dplyr)

# --- 2. Define Paths ---
base_dir <- normalizePath(file.path(here::here(), ".."))

objects_dir <- file.path(
  base_dir,
  "GitHub_Metabarcoding",
  "02_Analysis_Core",
  "objects"
)

subsets_dir <- file.path(
  base_dir,
  "GitHub_Metabarcoding",
  "02_Analysis_Core",
  "objects",
  "subsets"
)

if (!dir.exists(subsets_dir)) dir.create(subsets_dir, recursive = TRUE)

cat("\n======================================================\n")
cat("--- Integrating PR2 taxonomy and subsetting GreenAlgae ---\n")

# --- 3. Load Base Phyloseq Object and Taxonomy ---
cat("  > Loading biological-sample base phyloseq object...\n")

ps_base <- readRDS(
  file.path(objects_dir, "Algae_Phyloseq_Base.rds")
)

cat("  > Loading PR2 taxonomy object...\n")

tax_mat <- readRDS(
  file.path(objects_dir, "Algae_tax_final.rds")
)

cat("\nBase phyloseq object:\n")
print(ps_base)

cat("\nTaxonomy object dimensions:\n")
print(dim(tax_mat))

# --- 4. Clean and validate taxonomy matrix ---
cat("\n  > Cleaning taxonomy matrix...\n")

# PR2 / DADA2 expected levels:
expected_levels <- c(
  "Domain", "Supergroup", "Division", "Class",
  "Order", "Family", "Genus", "Species"
)

# Keep only first 8 columns if extra columns exist
if (ncol(tax_mat) >= 8) {
  tax_clean <- tax_mat[, 1:8, drop = FALSE]
  colnames(tax_clean) <- expected_levels
} else {
  stop("Taxonomy matrix has fewer than 8 columns. Check Algae_tax_final.rds.")
}

# Remove completely empty rows if any
tax_clean <- as.matrix(tax_clean)

# --- 5. Match taxonomy rows to phyloseq taxa ---
cat("  > Matching taxonomy rows to ASV names...\n")

asv_names <- taxa_names(ps_base)

# Case 1: taxonomy rownames already match ASV IDs
if (all(asv_names %in% rownames(tax_clean))) {
  
  tax_clean <- tax_clean[asv_names, , drop = FALSE]
  
} else if (nrow(tax_clean) == length(asv_names)) {
  
  # Case 2: taxonomy rows are in the same order but rownames are not ASV IDs
  warning(
    "Taxonomy rownames do not match phyloseq ASV IDs, ",
    "but the number of rows matches. Assigning phyloseq ASV IDs by order."
  )
  
  rownames(tax_clean) <- asv_names
  
} else {
  
  cat("\nFirst phyloseq taxa names:\n")
  print(head(asv_names, 10))
  
  cat("\nFirst taxonomy rownames:\n")
  print(head(rownames(tax_clean), 10))
  
  stop(
    "Taxonomy matrix cannot be matched to phyloseq taxa. ",
    "Check Script 04/05 ASV naming."
  )
}

if (!all(rownames(tax_clean) == taxa_names(ps_base))) {
  stop("Taxonomy rownames still do not match phyloseq taxa_names.")
}

# --- 6. Add taxonomy to phyloseq object ---
cat("  > Adding taxonomy to phyloseq object...\n")

tax_table(ps_base) <- tax_table(tax_clean)

ps_tax <- ps_base

cat("\nTaxonomized phyloseq object:\n")
print(ps_tax)

saveRDS(
  ps_tax,
  file.path(objects_dir, "Algae_Phyloseq_Master_Totalpool_TAX.rds")
)

# --- 7. Inspect taxonomy before target filtering ---
cat("\n  > Inspecting PR2 Division assignments...\n")

tax_df <- as.data.frame(as(tax_table(ps_tax), "matrix"))

cat("\nAvailable taxonomic ranks:\n")
print(colnames(tax_df))

cat("\nDivision summary:\n")
print(sort(table(tax_df$Division, useNA = "ifany"), decreasing = TRUE))

# --- 8. Filter target taxa ---
cat("\n  > Filtering target taxa: Chlorophyta:plas and Streptophyta:plas...\n")

target_divisions <- c(
  "Chlorophyta:plas",
  "Streptophyta:plas"
)

target_asvs <- rownames(tax_df)[
  !is.na(tax_df$Division) &
    tax_df$Division %in% target_divisions
]

if (length(target_asvs) == 0) {
  cat("\nNo ASVs matched the expected target divisions.\n")
  cat("Check the Division summary above and adjust target_divisions if needed.\n")
  stop("No target GreenAlgae ASVs retained.")
}

ps_algae <- prune_taxa(target_asvs, ps_tax)
ps_algae <- prune_taxa(taxa_sums(ps_algae) > 0, ps_algae)

cat("\nClean GreenAlgae phyloseq object:\n")
print(ps_algae)

cat("\n  > Total curated GreenAlgae / plastid ASVs:", ntaxa(ps_algae), "\n")
cat("  > Total samples:", nsamples(ps_algae), "\n")
cat("  > Total reads:", sum(sample_sums(ps_algae)), "\n")

# Save final clean master object
saveRDS(
  ps_algae,
  file.path(objects_dir, "Algae_Phyloseq_Master_Clean.rds")
)

# --- 9. Split by ecosystem using prune_samples ---
cat("\n  > Splitting data by Ecosystem...\n")

sample_meta <- data.frame(sample_data(ps_algae))

if (!"Ecosystem" %in% colnames(sample_meta)) {
  stop("Column 'Ecosystem' not found in sample_data(ps_algae).")
}

ecosystems <- unique(as.character(sample_meta$Ecosystem))
ecosystems <- ecosystems[!is.na(ecosystems)]

ps_ecosys_list <- list()

for (eco in ecosystems) {
  
  keep_samples <- rownames(sample_meta)[
    as.character(sample_meta$Ecosystem) == eco
  ]
  
  ps_sub <- prune_samples(keep_samples, ps_algae)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  ps_ecosys_list[[eco]] <- ps_sub
  
  cat("\n    -", eco, "\n")
  cat("      Samples:", nsamples(ps_sub), "\n")
  cat("      ASVs:", ntaxa(ps_sub), "\n")
  cat("      Reads:", sum(sample_sums(ps_sub)), "\n")
  
  # Save individual subset
  saveRDS(
    ps_sub,
    file.path(subsets_dir, paste0("Algae_", eco, ".rds"))
  )
}

# Save ecosystem list equivalent to CYA and Euka01
saveRDS(
  ps_ecosys_list,
  file.path(objects_dir, "Algae_Phyloseq_EcosysList_Totalpool_TAX.rds")
)

# --- 10. Final validation ---
cat("\n  > Final validation...\n")

cat("\nMaster clean object:\n")
cat("Samples:", nsamples(ps_algae), "\n")
cat("ASVs:", ntaxa(ps_algae), "\n")
cat("First sample names:\n")
print(head(sample_names(ps_algae), 20))

cat("\nEcosystem table:\n")
print(table(sample_data(ps_algae)$Ecosystem))

cat("\n[SUCCESS] Script 06 completed.\n")
cat("Saved:\n")
cat(" - Algae_Phyloseq_Master_Totalpool_TAX.rds\n")
cat(" - Algae_Phyloseq_Master_Clean.rds\n")
cat(" - Algae_Phyloseq_EcosysList_Totalpool_TAX.rds\n")
cat(" - Ecosystem subsets in objects/subsets/\n")
cat("======================================================\n")