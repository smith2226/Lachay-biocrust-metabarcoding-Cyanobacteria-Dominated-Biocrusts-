# ==============================================================================
# Script 03: Full Contamination Filtering (Cyanobacteria 16S) - MATRIX FIXED
# Project: Biological Soil Crust Peru
# Purpose: Sequential subtraction of blanks & 2/4 Rule (Bulletproof Version)
# ==============================================================================

# --- 1. Load Libraries ---
library(dplyr)
library(tidyr)
library(here)
# --- 2. Define Paths ---
input_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
meta_path  <- here("Zenodo_Metabarcoding", "01_Data", "Metadata", "Pierr_metadata.csv")
diag_path  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Cyano_Filtering")

if(!dir.exists(diag_path)) dir.create(diag_path, recursive = TRUE)

# --- 3. Load Metadata ---
if (!file.exists(meta_path)) stop("Metadata file not found! Check path.")
exdata <- read.csv(meta_path, stringsAsFactors = FALSE)

# --- 4. Helper Functions ---
get_stats <- function(otu_matrix) {
  c(Reads = sum(otu_matrix, na.rm = TRUE), OTUs = sum(colSums(otu_matrix, na.rm = TRUE) > 0))
}

rep.groups2 <- function(x){
  r2 <- apply(x, 2, function(col_otu) {
    # Safe sum ignoring any NAs just in case
    if(sum(col_otu > 0, na.rm = TRUE) < 2) rep(0, length(col_otu)) else col_otu
  })
  if (is.vector(r2)) r2 <- matrix(r2, nrow=nrow(x), ncol=ncol(x))
  dimnames(r2) <- dimnames(x)
  r2[is.na(r2)] <- 0 # Ensure absolute zero NAs
  return(as.matrix(r2))
}

# --- 5. Main Processing Loop ---
pooling_types <- c("independent", "totalpool", "pseudopool")
qc_report_list <- list()

for (type in pooling_types) {
  cat("\n======================================================\n")
  cat("--- Running Decontamination for Strategy:", toupper(type), "---\n")
  
  # A. Load LULU Object
  input_file <- file.path(input_dir, paste0("Cyano_LULU_curated_", type, ".rds"))
  if (!file.exists(input_file)) {
    warning(paste("File not found, skipping:", input_file))
    next
  }
  
  lulu_obj <- readRDS(input_file)
  otu_table <- as.matrix(t(lulu_obj$curated_table))
  rownames(otu_table) <- sub("_F_filt\\.fastq\\.gz", "", rownames(otu_table))
  stats_initial <- get_stats(otu_table)
  
  # Map Metadata directly to the rows of otu_table
  rep_bases <- sub("_rep[1-4]$", "", rownames(otu_table))
  
  rep_type <- exdata$TypeSample_Treatment[match(rep_bases, exdata$SampleID)]
  rep_type[is.na(rep_type)] <- "Unknown"
  
  rep_ext_batch <- exdata$IDseq_ExtractBlank[match(rep_bases, exdata$SampleID)]
  rep_ext_batch[is.na(rep_ext_batch) | rep_ext_batch == ""] <- "No_Batch"
  
  rep_fld_batch <- exdata$IDseq_FieldControl[match(rep_bases, exdata$SampleID)]
  rep_fld_batch[is.na(rep_fld_batch) | rep_fld_batch == ""] <- "No_Batch"
  
  rep_is_bio <- exdata$BiologicalSample[match(rep_bases, exdata$SampleID)]
  
  # B. STEP 1: NTC Subtraction
  cat("  > Step 1: NTC Subtraction...\n")
  ntc_reps <- rownames(otu_table)[rep_type == "Negative_LAB"]
  if(length(ntc_reps) > 0) {
    ntc_max <- apply(otu_table[ntc_reps, , drop = FALSE], 2, max, na.rm = TRUE)
    ntc_max[is.infinite(ntc_max)] <- 0
    otu_s1 <- sweep(otu_table, 2, ntc_max, FUN = "-")
    otu_s1 <- pmax(otu_s1, 0) # Convert negatives to 0 securely
  } else {
    otu_s1 <- otu_table
  }
  stats_s1 <- get_stats(otu_s1)
  
  # C. STEP 2: Extraction Blank Subtraction (By Batch)
  cat("  > Step 2: Extraction Blank Subtraction...\n")
  otu_s2 <- otu_s1
  unique_ext_batches <- unique(rep_ext_batch[rep_ext_batch != "No_Batch"])
  
  for(batch in unique_ext_batches) {
    batch_blanks <- rownames(otu_s2)[rep_ext_batch == batch & rep_type == "Extraction_Blank"]
    if(length(batch_blanks) > 0) {
      blank_max <- apply(otu_s2[batch_blanks, , drop = FALSE], 2, max, na.rm = TRUE)
      blank_max[is.infinite(blank_max)] <- 0
      
      batch_targets <- rownames(otu_s2)[rep_ext_batch == batch]
      subtracted <- sweep(otu_s2[batch_targets, , drop = FALSE], 2, blank_max, FUN = "-")
      otu_s2[batch_targets, ] <- pmax(subtracted, 0)
    }
  }
  stats_s2 <- get_stats(otu_s2)
  
  # D. STEP 3: Field Blank Subtraction (By Batch)
  cat("  > Step 3: Field Blank Subtraction...\n")
  otu_s3 <- otu_s2
  unique_fld_batches <- unique(rep_fld_batch[rep_fld_batch != "No_Batch"])
  
  for(batch in unique_fld_batches) {
    batch_blanks <- rownames(otu_s3)[rep_fld_batch == batch & rep_type == "Field_Blank"]
    if(length(batch_blanks) > 0) {
      blank_max <- apply(otu_s3[batch_blanks, , drop = FALSE], 2, max, na.rm = TRUE)
      blank_max[is.infinite(blank_max)] <- 0
      
      batch_targets <- rownames(otu_s3)[rep_fld_batch == batch]
      subtracted <- sweep(otu_s3[batch_targets, , drop = FALSE], 2, blank_max, FUN = "-")
      otu_s3[batch_targets, ] <- pmax(subtracted, 0)
    }
  }
  stats_s3 <- get_stats(otu_s3)
  
  # E. STEP 4: 2/4 Rule
  cat("  > Step 4: Biological Filter & Replicate Consistency (2/4)...\n")
  bio_reps <- rownames(otu_s3)[rep_is_bio == "Yes"]
  otu_bio <- otu_s3[bio_reps, , drop = FALSE]
  
  sample_bases <- sub("_rep[1-4]$", "", rownames(otu_bio))
  list_reps <- split(as.data.frame(otu_bio), sample_bases)
  
  otu_s4_list <- lapply(list_reps, rep.groups2)
  otu_final <- do.call(rbind, otu_s4_list)
  rownames(otu_final) <- sub("^.*\\.", "", rownames(otu_final)) # Clean names
  
  # Clean empty OTUs
  otu_final <- otu_final[, colSums(otu_final, na.rm = TRUE) > 0, drop = FALSE]
  stats_s4 <- get_stats(otu_final)
  
  # F. Export Master Table
  saveRDS(otu_final, file.path(output_dir, paste0("Cyano_final_master_", type, ".rds")))
  
  cat(sprintf("  [RESULTS %s]\n", toupper(type)))
  cat(sprintf("   - Final Bio:    %d Reads | %d OTUs\n", stats_s4[1], stats_s4[2]))
  
  qc_report_list[[type]] <- data.frame(
    Pooling_Strategy = type,
    Step = c("01_Initial", "02_After_NTC", "03_After_Ext_Blank", "04_After_Field_Blank", "05_Final_Bio_2of4_Rule"),
    Total_Reads = c(stats_initial[1], stats_s1[1], stats_s2[1], stats_s3[1], stats_s4[1]),
    Total_OTUs = c(stats_initial[2], stats_s1[2], stats_s2[2], stats_s3[2], stats_s4[2])
  )
}

# --- 6. Master QC Report ---
qc_report_df <- bind_rows(qc_report_list)
write.csv(qc_report_df, file.path(diag_path, "Cyano_Contamination_Filtering_QC_Report.csv"), row.names = FALSE)

cat("\n======================================================\n")
cat("--- Script 03 Complete! Corrected .rds files generated. ---\n")