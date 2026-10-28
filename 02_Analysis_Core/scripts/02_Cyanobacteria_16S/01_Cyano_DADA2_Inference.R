# ==============================================================================
# Script 01: DADA2 ASV Inference (Cyanobacteria 16S) - EcoPC Version
# Project: Biological Soil Crust Peru
# Purpose: From trimmed fastq files to Chimera-free ASV tables (NovaSeq Corrected)
# ==============================================================================

#### Moving into R for DADA2 processing...
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(dada2); packageVersion("dada2") 
library(ShortRead); packageVersion("ShortRead") 
library(dplyr); packageVersion("dplyr") 
library(tidyr); packageVersion("tidyr") 
library(Hmisc); packageVersion("Hmisc") 
library(ggplot2); packageVersion("ggplot2") 
library(plotly); packageVersion("plotly") 
library(here)

# --- 2. Paths (Windows format for RStudio) ---
path <- here("Zenodo_Metabarcoding", "01_Data", "Raw", "Cyano_Green_Mixed", "cyanoseqfiles", "Cyanotrimmed")

if (!dir.exists(path)) stop("Path not found! Check your Zenodo_Metabarcoding directory structure.")

setwd(path)

output_path <- here("GitHub_Metabarcoding", "02_Analysis_Core", "objects")
diag_path   <- here("GitHub_Metabarcoding", "02_Analysis_Core", "results", "quality_reports")

if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if(!dir.exists(diag_path)) dir.create(diag_path, recursive = TRUE)

# --- 3. File List & Sample Extraction ---
# Matching your actual file structure: S001.rep1.trim1.fq.gz
fnFs <- sort(list.files(path, pattern = "\\.trim1\\.fq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "\\.trim2\\.fq\\.gz$", full.names = TRUE))

# Verify equal number of Forward and Reverse files
stopifnot(length(fnFs) == length(fnRs))

# Extraction logic for S001.rep1 format -> S001_rep1
sample.names <- sapply(basename(fnFs), function(x) {
  matches <- regmatches(x, regexec("^(S[0-9]+)\\.rep([0-9]+)", x))
  if (length(matches[[1]]) >= 3) {
    paste0(matches[[1]][2], "_rep", matches[[1]][3])
  } else { NA }
})

# --- 4. Filtering and Trimming ---
filt_path <- file.path(path, "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## --- 4.1 Checking and Saving Quality Reports ---
pdf(file.path(diag_path, "Cyano_Initial_Quality_Full_Report.pdf"), width = 12, height = 15)
for(i in seq(1, length(fnFs), by = 10)) {
  end <- min(i + 9, length(fnFs))
  print(plotQualityProfile(fnFs[i:end]) + ggtitle(paste("Forward Reads - Samples", i, "to", end)))
  print(plotQualityProfile(fnRs[i:end]) + ggtitle(paste("Reverse Reads - Samples", i, "to", end)))
}
dev.off()

# Filtering: truncLen is set to 0 as reads are already short or as standard
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = FALSE)

# Remove samples with zero reads after filtering
keep <- out[,"reads.out"] > 0
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]
sample.names <- sample.names[keep]; out <- out[keep,]

# ==============================================================================
# 5. Error Learning: Default vs. Hacked
# ==============================================================================

# 5.1 DEFAULT ERROR 
errF_default <- learnErrors(filtFs, multithread = TRUE, MAX_CONSIST = 20)
errR_default <- learnErrors(filtRs, multithread = TRUE, MAX_CONSIST = 20)

pdf(file.path(diag_path, "Cyano_Error_Plots_DEFAULT.pdf"))
print(plotErrors(errF_default, nominalQ = TRUE) + ggtitle("Forward Errors - DEFAULT"))
print(plotErrors(errR_default, nominalQ = TRUE) + ggtitle("Reverse Errors - DEFAULT"))
dev.off()

# --- 5.2 NovaSeq Error Model (JacobRPrice & Guillem Salazar Solution) ---
library(magrittr)
library(dplyr)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } 
    } 
  } 
  
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  return(err)
}

errF <- learnErrors(filtFs, multithread = TRUE, errorEstimationFunction = loessErrfun_mod, verbose = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE, errorEstimationFunction = loessErrfun_mod, verbose = TRUE)

pdf(file.path(diag_path, "Cyano_Error_Plots_HACKED.pdf"))
print(plotErrors(errF, nominalQ = TRUE) + ggtitle("Forward Errors - NOVASEQ HACK"))
print(plotErrors(errR, nominalQ = TRUE) + ggtitle("Reverse Errors - NOVASEQ HACK"))
dev.off()

# ==============================================================================
# 6. DADA2 Core: NP, PP, and Total Pooling
# ==============================================================================
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Running Independent, Pseudo-pooling, and Total Pooling strategies
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFpsPPs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

dadaFPPs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

save.image(file.path(output_path, "Cyano_dada2_analysis.RData"))

# ==============================================================================
# 7. Merging & Chimera Removal
# ==============================================================================

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
mergersPP <- mergePairs(dadaFPPs, derepFs, dadaRPPs, derepRs, verbose=TRUE)
mergers_psPP <- mergePairs(dadaFpsPPs, derepFs, dadaRpsPPs, derepRs, verbose=TRUE)

## Making sequence table
P1seqtab <- makeSequenceTable(mergers)
P1seqtabPP <- makeSequenceTable(mergersPP)
P1seqtabpsPP <- makeSequenceTable(mergers_psPP)

## Removing Chimeras
P1seqtab.nochim <- removeBimeraDenovo(P1seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
P1seqtabPP.nochim <- removeBimeraDenovo(P1seqtabPP, method="consensus", multithread=TRUE, verbose=TRUE)
P1seqtabpsPP.nochim <- removeBimeraDenovo(P1seqtabpsPP, method="consensus", multithread=TRUE, verbose=TRUE)

save.image(file.path(output_path, "Cyano_dada2_analysis.RData"))

# ==============================================================================
# 8. Tracking Pipeline Progress (QC)
# ==============================================================================

getN <- function(x) sum(getUniques(x))

# 8.1 Tracking - Independent (No pooling)
track_NP <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
                  sapply(mergers, getN), rowSums(P1seqtab.nochim))
colnames(track_NP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_NP) <- sample.names

# 8.2 Tracking - Total Pooling
track_PP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), 
                  sapply(mergersPP, getN), rowSums(P1seqtabPP.nochim))
colnames(track_PP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_PP) <- sample.names

# 8.3 Tracking - Pseudo-pooling
track_psPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), 
                    sapply(mergers_psPP, getN), rowSums(P1seqtabpsPP.nochim))
colnames(track_psPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_psPP) <- sample.names

# Exporting all tracking tables
write.csv(track_NP, file.path(diag_path, "Cyano_Tracking_Independent.csv"))
write.csv(track_PP, file.path(diag_path, "Cyano_Tracking_TotalPooling.csv"))
write.csv(track_psPP, file.path(diag_path, "Cyano_Tracking_PseudoPooling.csv"))

# ==============================================================================
# 9. Final Exports
# ==============================================================================

saveRDS(P1seqtab.nochim, file.path(output_path, "Cyano_seqtab_nochim_independent.rds"))
saveRDS(P1seqtabPP.nochim, file.path(output_path, "Cyano_seqtab_nochim_totalpool.rds"))
saveRDS(P1seqtabpsPP.nochim, file.path(output_path, "Cyano_seqtab_nochim_pseudopool.rds"))

save.image(file.path(output_path, "Cyano_DADA2_final_complete_workspace.RData"))

# Exporting summary table
asv_summary <- data.frame(
  Method = c("Independent", "Total_Pooling", "Pseudo_Pooling"),
  Total_ASVs = c(ncol(P1seqtab.nochim), ncol(P1seqtabPP.nochim), ncol(P1seqtabpsPP.nochim)),
  Total_Reads = c(sum(P1seqtab.nochim), sum(P1seqtabPP.nochim), sum(P1seqtabpsPP.nochim))
)

write.csv(asv_summary, file.path(diag_path, "Cyano_ASV_Method_Comparison_Summary.csv"), row.names = FALSE)

cat("\nStep 01 finished successfully. All objects saved to:", output_path, "\n")