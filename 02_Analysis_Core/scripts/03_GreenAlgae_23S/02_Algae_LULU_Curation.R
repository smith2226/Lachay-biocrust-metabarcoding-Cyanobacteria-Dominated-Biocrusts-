# ==============================================================================
# Script 02: LULU ASV Curation (Green Algae) - ALL POOLING STRATEGIES
# Project: Biological Soil Crust Peru
# Purpose: Collapse redundant/erroneous ASVs for Independent, Totalpool, and Pseudopool
# ==============================================================================

# --- 1. Load Libraries ---
library(dada2)
library(lulu)
library(dplyr)
library(Biostrings)
library(here)
# --- 2. Define Paths ---
input_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
blast_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Algae_LULU")

if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
if(!dir.exists(blast_dir)) dir.create(blast_dir, recursive=TRUE)
setwd(blast_dir)
cat("\n======================================================\n")
cat("--- Running LULU Curation for Green Algae (ALL METHODS) ---\n")

# --- 3. Define Strategies ---
strategies <- c("independent", "totalpool", "pseudopool")

# Dataframe to store summary results
lulu_summary <- data.frame(
  Strategy = character(),
  Original_ASVs = numeric(),
  Curated_ASVs = numeric(),
  Discarded_ASVs = numeric(),
  stringsAsFactors = FALSE
)

# --- 4. Loop Through Each Strategy ---
for (strat in strategies) {
  cat(paste0("\n  > [", toupper(strat), "] Starting processing...\n"))
  
  # A. Load Sequence Table
  seqtab_file <- file.path(input_dir, paste0("Algae_seqtab_nochim_", strat, ".rds"))
  if (!file.exists(seqtab_file)) {
    cat(paste0("    [!] File not found: ", seqtab_file, " - Skipping.\n"))
    next
  }
  seqtab <- readRDS(seqtab_file)
  
  # B. Format for LULU
  asv_table <- as.data.frame(t(seqtab))
  seqs <- rownames(asv_table)
  asv_ids <- paste0("ASV_", 1:nrow(asv_table))
  rownames(asv_table) <- asv_ids
  
  # C. Create FASTA
  fasta_file <- file.path(blast_dir, paste0("Algae_ASVs_", strat, ".fasta"))
  uniquesToFasta(seqtab, fasta_file, ids = asv_ids)
  
  # D. Run BLAST
  setwd(blast_dir)
  cat(paste0("    - Running BLASTn for ", strat, "...\n"))
  
  # Build Database
  cmd_makedb <- paste0("makeblastdb -in Algae_ASVs_", strat, ".fasta -parse_seqids -dbtype nucl")
  system(cmd_makedb, ignore.stdout = TRUE)
  
  # Run Alignment
  match_file <- paste0("Algae_match_list_", strat, ".txt")
  cmd_blast <- paste0("blastn -db Algae_ASVs_", strat, ".fasta -outfmt \"6 qseqid sseqid pident\" -out ", match_file, " -qcov_hsp_perc 80 -perc_identity 84 -query Algae_ASVs_", strat, ".fasta")
  system(cmd_blast)
  
  # E. Run LULU
  cat(paste0("    - Running LULU algorithm for ", strat, "...\n"))
  match_list <- read.table(match_file, header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
  curated_result <- lulu(asv_table, match_list)
  
  # F. Reconstruct and Save
  curated_table <- curated_result$curated_table
  names(seqs) <- asv_ids
  curated_seqs <- seqs[rownames(curated_table)]
  
  seqtab_curated <- t(curated_table)
  colnames(seqtab_curated) <- curated_seqs
  
  saveRDS(seqtab_curated, file.path(output_dir, paste0("Algae_seqtab_nochim_", strat, "_LULU.rds")))
  saveRDS(curated_result, file.path(blast_dir, paste0("Algae_LULU_Result_", strat, ".rds")))
  
  # G. Record Stats
  lulu_summary <- rbind(lulu_summary, data.frame(
    Strategy = strat,
    Original_ASVs = nrow(asv_table),
    Curated_ASVs = curated_result$curated_count,
    Discarded_ASVs = curated_result$discarded_count
  ))
  
  cat(paste0("    [OK] ", toupper(strat), " curated: ", curated_result$curated_count, " ASVs retained.\n"))
}

# --- 5. Export Final Summary ---
write.csv(lulu_summary, file.path(blast_dir, "Algae_LULU_Method_Comparison_Summary.csv"), row.names = FALSE)

cat("\n======================================================\n")
cat("[SUCCESS] LULU Curation Completed for ALL Methods!\n")
print(lulu_summary)
cat("======================================================\n")