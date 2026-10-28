# ==============================================================================
# Script 03: Full Contamination Filtering (Eukaryotes 18S)
# Project: Biological Soil Crust Peru
# Purpose: Sequential subtraction of blanks & QC reporting (All Strategies)
# ==============================================================================

# --- 1. Load Libraries ---
library(dplyr)
library(tidyr)
library(here)
# --- 2. Define Paths (Zenodo FAIR Structure) ---
input_dir  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
output_dir <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
meta_path  <- here("Zenodo_Metabarcoding", "01_Data", "Metadata", "Pierr_metadata.csv")
lulu_dir   <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Euka_LULU")
diag_path  <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports", "Euka_Filtering")

if(!dir.exists(diag_path)) dir.create(diag_path, recursive = TRUE)

# --- 3. Load Metadata ---
if (!file.exists(meta_path)) stop("Metadata file not found! Check path.")
exdata <- read.csv(meta_path, stringsAsFactors = FALSE)

# --- 4. Filtering Functions ---

ntc.change <- function(otu_table, ntc_replicate_ids){
  ntc_rows <- intersect(rownames(otu_table), ntc_replicate_ids)
  if (length(ntc_rows) == 0) return(otu_table)
  max_ntc <- apply(otu_table[ntc_rows, , drop = FALSE], 2, function(y) max(y, na.rm = TRUE))
  max_ntc[is.infinite(max_ntc)] <- 0
  otu_cleaned <- sweep(otu_table, 2, max_ntc, FUN = "-")
  return(pmax(otu_cleaned, 0))
}

blank.change <- function(df_batch, otu_cols, type_col, label) {
  blank_rows <- rownames(df_batch)[which(df_batch[[type_col]] == label)]
  otu_data <- as.matrix(df_batch[, otu_cols, drop = FALSE])
  mode(otu_data) <- "numeric"
  if (length(blank_rows) == 0) return(otu_data)
  max_blank <- apply(otu_data[blank_rows, , drop = FALSE], 2, function(y) max(y, na.rm = TRUE))
  max_blank[is.infinite(max_blank)] <- 0
  cleaned <- sweep(otu_data, 2, max_blank, FUN = "-")
  return(pmax(cleaned, 0))
}

rep.groups2 <- function(x){
  r2 <- apply(x, 2, function(col_otu) {
    if(sum(col_otu > 0) < 2) rep(0, length(col_otu)) else col_otu
  })
  if (is.vector(r2)) r2 <- matrix(r2, nrow=nrow(x), ncol=ncol(x))
  
  # PARCHE 1: Forzar conservación de nombres de filas y columnas
  dimnames(r2) <- dimnames(x)
  return(as.matrix(r2))
}

get_stats <- function(otu_table) {
  c(Reads = sum(otu_table), OTUs = sum(colSums(otu_table) > 0))
}

# --- 5. Main Processing Loop & QC Tracking ---
pooling_types <- c("independent", "totalpool", "pseudopool")
qc_report_list <- list()

for (type in pooling_types) {
  cat("\n======================================================\n")
  cat("--- Running Decontamination for Strategy:", toupper(type), "---\n")
  
  # A. Load Table & Index
  input_file <- file.path(input_dir, paste0("Euka_seqtab_curated_", type, ".rds"))
  index_file <- file.path(lulu_dir, paste0("Euka_indices_curated_", type, ".csv"))
  
  if (!file.exists(input_file)) {
    warning(paste("File not found, skipping:", input_file))
    next
  }
  
  otu_table <- readRDS(input_file)
  index_data <- read.csv(index_file, stringsAsFactors = FALSE)
  otu_cols <- colnames(otu_table)
  stats_initial <- get_stats(otu_table)
  
  # B. STEP 1: NTC 
  cat("  > Step 1: NTC Subtraction...\n")
  ntc_samples <- exdata$SampleID[which(exdata$TypeSample_Treatment == "Negative_LAB")]
  ntc_ids <- unique(unlist(lapply(ntc_samples, function(s) paste0(s, "_rep", 1:4))))
  otu_s1 <- ntc.change(otu_table, ntc_ids)
  stats_s1 <- get_stats(otu_s1)
  
  # C. STEP 2: Extraction Blank 
  cat("  > Step 2: Extraction Blank Subtraction...\n")
  z_merged <- merge(otu_s1, index_data[, c("full", "sample")], by.x = "row.names", by.y = "full")
  z_final  <- merge(z_merged, exdata[, c("SampleID", "IDseq_ExtractBlank", "IDseq_FieldControl", "TypeSample_Treatment", "BiologicalSample")], by.x = "sample", by.y = "SampleID")
  rownames(z_final) <- z_final$Row.names
  
  z_final$IDseq_ExtractBlank[is.na(z_final$IDseq_ExtractBlank) | z_final$IDseq_ExtractBlank == ""] <- "No_Batch"
  
  list_eb <- split(z_final, z_final$IDseq_ExtractBlank)
  otu_eb_list <- lapply(list_eb, blank.change, otu_cols = otu_cols, type_col = "TypeSample_Treatment", label = "Extraction_Blank")
  
  otu_s2 <- do.call(rbind, otu_eb_list)
  # PARCHE 2: Limpiar los nombres de filas que R ensucia al hacer rbind
  rownames(otu_s2) <- sub("^.*\\.", "", rownames(otu_s2))
  otu_s2 <- otu_s2[rownames(otu_table), ] 
  stats_s2 <- get_stats(otu_s2)
  
  # D. STEP 3: Field Blank 
  cat("  > Step 3: Field Blank Subtraction...\n")
  z_final_s3 <- z_final
  z_final_s3[, otu_cols] <- otu_s2[rownames(z_final_s3), otu_cols] 
  
  z_final_s3$IDseq_FieldControl[is.na(z_final_s3$IDseq_FieldControl) | z_final_s3$IDseq_FieldControl == ""] <- "No_Batch"
  
  list_fb <- split(z_final_s3, z_final_s3$IDseq_FieldControl)
  otu_fb_list <- lapply(list_fb, blank.change, otu_cols = otu_cols, type_col = "TypeSample_Treatment", label = "Field_Blank")
  
  otu_s3 <- do.call(rbind, otu_fb_list)
  rownames(otu_s3) <- sub("^.*\\.", "", rownames(otu_s3))
  otu_s3 <- otu_s3[rownames(otu_table), ] 
  stats_s3 <- get_stats(otu_s3)
  
  # E. STEP 4: 2/4 Rule
  cat("  > Step 4: Biological Filter & Replicate Consistency (2/4)...\n")
  bio_ids <- z_final$Row.names[which(z_final$BiologicalSample == "Yes")]
  otu_bio <- otu_s3[intersect(rownames(otu_s3), bio_ids), ]
  
  sample_bases <- sub("_rep[1-4]$", "", rownames(otu_bio))
  list_reps <- split(as.data.frame(otu_bio), sample_bases)
  
  otu_s4_list <- lapply(list_reps, rep.groups2)
  otu_final <- do.call(rbind, otu_s4_list)
  rownames(otu_final) <- sub("^.*\\.", "", rownames(otu_final))
  
  # Eliminar OTUs con 0 abundancia total
  otu_final <- otu_final[, colSums(otu_final) > 0, drop = FALSE]
  stats_s4 <- get_stats(otu_final)
  
  # F. Export Master Table
  saveRDS(otu_final, file.path(output_dir, paste0("Euka_final_master_", type, ".rds")))
  
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
write.csv(qc_report_df, file.path(diag_path, "Euka_Contamination_Filtering_QC_Report.csv"), row.names = FALSE)

save.image(file.path(output_dir, "Euka_Script03_Filtering_final.RData"))
cat("\n======================================================\n")
cat("--- Script 03 Complete! Corrected .rds files generated. ---\n")
