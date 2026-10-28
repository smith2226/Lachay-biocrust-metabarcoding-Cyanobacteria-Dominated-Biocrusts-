# ==============================================================================
# Script 01: Baseline Metabarcoding Summary for Paper01
# Paper: Lomas Biomass & Nutrient Addition
# Purpose:
#   Generate final audit, taxonomic summaries, alpha diversity,
#   supplementary figures, and text blocks for the initial T3 inoculum.
# ==============================================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(vegan)

# --- 1. Paths (relative to working dir: 03_Paper01_Lomas_Biomass/) ---
DATA_DIR <- "data"
tab_dir  <- file.path("results", "tables", "Metabarcoding_Baseline")
fig_dir  <- file.path("results", "figures", "Metabarcoding_Baseline")
text_dir <- file.path("results", "text")

dir.create(tab_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)

# --- PR2 label cleaning (consistent with P02/P03/P04) ---
clean_pr2_label <- function(x) {
  x <- gsub(":plas$", "", x)
  x <- gsub("_(XX|X)$", "", x)
  x <- gsub("_[A-Z]{2,}$", "", x)
  x <- gsub("_", " ", x)
  x <- str_trim(x)
  x[is.na(x) | x == ""] <- "Unclassified"
  x
}

# --- 2. Load final Paper01 objects ---
files <- list(
  CYA = file.path(DATA_DIR, "CYA_T3_Lachay_baseline_final.rds"),
  GreenAlgae = file.path(DATA_DIR, "GreenAlgae_T3_Lachay_baseline_final.rds"),
  Euka01 = file.path(DATA_DIR, "Euka01_T3_Lachay_baseline_final.rds")
)

missing_files <- files[!file.exists(unlist(files))]

if (length(missing_files) > 0) {
  print(missing_files)
  stop("Some baseline RDS files are missing. Run marker distribution/audit first.")
}

ps_cya <- readRDS(files$CYA)
ps_green <- readRDS(files$GreenAlgae)
ps_euka <- readRDS(files$Euka01)

# --- 3. Validate samples ---
target_ids <- c("S082", "S083", "S084", "S085", "S086")

validate_ps <- function(ps, marker) {
  
  if (nsamples(ps) != 5) {
    stop(marker, " should contain 5 biological samples, but has ", nsamples(ps))
  }
  
  missing_ids <- setdiff(target_ids, sample_names(ps))
  
  if (length(missing_ids) > 0) {
    stop(marker, " is missing samples: ", paste(missing_ids, collapse = ", "))
  }
  
  cat(marker, "validated:", nsamples(ps), "samples,", ntaxa(ps), "ASVs\n")
}

validate_ps(ps_cya, "CYA")
validate_ps(ps_green, "GreenAlgae")
validate_ps(ps_euka, "Euka01")

# --- 4. Audit summary ---
audit_phyloseq <- function(ps, marker, object_name) {
  
  ss <- sample_sums(ps)
  
  data.frame(
    Marker = marker,
    Object_used = object_name,
    Samples = nsamples(ps),
    Biological_samples = paste(sample_names(ps), collapse = "; "),
    ASVs = ntaxa(ps),
    Total_reads = sum(ss),
    Min_reads_per_sample = min(ss),
    Median_reads_per_sample = median(ss),
    Max_reads_per_sample = max(ss),
    stringsAsFactors = FALSE
  )
}

audit_table <- bind_rows(
  audit_phyloseq(ps_cya, "CYA_16S_Cyanobacteria", "Cyano_Phyloseq_EcosysList_Totalpool_TAX.rds"),
  audit_phyloseq(ps_green, "GreenAlgae_G500_G800", "Algae_Phyloseq_EcosysList_Totalpool_TAX.rds"),
  audit_phyloseq(ps_euka, "Euka01_18S_Eukaryotes", "Euka_Phyloseq_EcosysList_Totalpool_TAX.rds")
)

write_csv(
  audit_table,
  file.path(tab_dir, "S9_Metabarcoding_Audit_FinalObjects.csv")
)

print(audit_table)

# --- 5. Taxonomic summaries ---
make_tax_summary <- function(ps, marker, rank_name = "Genus") {
  
  ranks <- rank_names(ps)
  
  if (!(rank_name %in% ranks)) {
    warning("Rank ", rank_name, " not found in ", marker, ". Using last available rank.")
    rank_name <- tail(ranks, 1)
  }
  
  otu <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) otu <- as.data.frame(t(otu))
  otu$ASV <- rownames(otu)
  
  tax <- as.data.frame(as(tax_table(ps), "matrix"))
  tax$ASV <- rownames(tax)
  
  otu_long <- otu %>%
    pivot_longer(
      cols = -ASV,
      names_to = "Sample",
      values_to = "Reads"
    ) %>%
    group_by(ASV) %>%
    summarise(
      Total_reads = sum(Reads, na.rm = TRUE),
      Detection_frequency = sum(Reads > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  out <- otu_long %>%
    left_join(tax, by = "ASV") %>%
    mutate(
      Marker = marker,
      Rank_used = rank_name,
      Taxon = clean_pr2_label(as.character(.data[[rank_name]]))
    ) %>%
    group_by(Marker, Rank_used, Taxon) %>%
    summarise(
      ASVs = n(),
      Total_reads = sum(Total_reads, na.rm = TRUE),
      Mean_detection_frequency = mean(Detection_frequency),
      .groups = "drop"
    ) %>%
    group_by(Marker) %>%
    mutate(
      Relative_abundance_percent = Total_reads / sum(Total_reads) * 100
    ) %>%
    arrange(Marker, desc(Relative_abundance_percent))
  
  return(out)
}

tax_cya <- make_tax_summary(ps_cya, "CYA_16S_Cyanobacteria", rank_name = "Genus")
tax_green <- make_tax_summary(ps_green, "GreenAlgae_G500_G800", rank_name = "Genus")
tax_euka <- make_tax_summary(ps_euka, "Euka01_18S_Eukaryotes", rank_name = "Class")

write_csv(tax_cya, file.path(tab_dir, "S5_CYA_TaxonomicSummary_Genus.csv"))
write_csv(tax_green, file.path(tab_dir, "S6_GreenAlgae_TaxonomicSummary_Genus.csv"))
write_csv(tax_euka, file.path(tab_dir, "S7_Euka01_TaxonomicSummary_Class.csv"))

# --- 6. Alpha diversity ---
alpha_summary <- function(ps, marker) {
  
  otu <- as.data.frame(otu_table(ps))
  if (taxa_are_rows(ps)) otu <- as.data.frame(t(otu))
  
  reads <- rowSums(otu)
  observed <- vegan::specnumber(otu)
  shannon <- vegan::diversity(otu, index = "shannon")
  simpson <- vegan::diversity(otu, index = "simpson")
  pielou <- ifelse(observed > 1, shannon / log(observed), NA)
  
  data.frame(
    Marker = marker,
    Sample = rownames(otu),
    Reads = reads,
    Observed_ASVs = observed,
    Shannon = shannon,
    Simpson = simpson,
    Pielou_evenness = pielou,
    stringsAsFactors = FALSE
  )
}

alpha_table <- bind_rows(
  alpha_summary(ps_cya, "CYA_16S_Cyanobacteria"),
  alpha_summary(ps_green, "GreenAlgae_G500_G800"),
  alpha_summary(ps_euka, "Euka01_18S_Eukaryotes")
)

write_csv(
  alpha_table,
  file.path(tab_dir, "S8_AlphaDiversity_T3_Baseline.csv")
)

# --- 7. Supplementary barplot function ---
plot_composition <- function(ps, marker_title, rank_name = "Genus", top_n = 12,
                             palette_base = "Paired") {
  
  if (!(rank_name %in% rank_names(ps))) {
    rank_name <- tail(rank_names(ps), 1)
  }
  
  pd <- psmelt(ps)
  
  pd[[rank_name]] <- clean_pr2_label(as.character(pd[[rank_name]]))
  
  pd_agg <- pd %>%
    group_by(Sample, .data[[rank_name]]) %>%
    summarise(
      Abundance = sum(Abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Sample) %>%
    mutate(
      RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)
    ) %>%
    ungroup()
  
  top_taxa <- pd_agg %>%
    group_by(.data[[rank_name]]) %>%
    summarise(
      MeanRel = mean(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(MeanRel)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[rank_name]])
  
  pd_final <- pd_agg %>%
    mutate(
      Taxon_plot = ifelse(.data[[rank_name]] %in% top_taxa,
                          .data[[rank_name]], "Other")
    ) %>%
    group_by(Sample, Taxon_plot) %>%
    summarise(
      RelAbundance = sum(RelAbundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  pd_final$Taxon_plot <- factor(pd_final$Taxon_plot, levels = c(top_taxa, "Other"))
  
  nlev <- length(levels(pd_final$Taxon_plot))
  brewer_n <- max(3, min(nlev, 12))
  pal <- colorRampPalette(brewer.pal(brewer_n, palette_base))(nlev)
  pal[nlev] <- "#D3D3D3"
  
  p <- ggplot(pd_final, aes(x = Sample, y = RelAbundance, fill = Taxon_plot)) +
    geom_col(color = "black", linewidth = 0.15) +
    scale_y_continuous(
      labels = scales::percent_format(),
      expand = c(0, 0.005)
    ) +
    scale_fill_manual(values = pal, name = rank_name) +
    labs(
      title = marker_title,
      x = "Initial biocrust samples",
      y = "Relative abundance (%)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8, face = "italic"),
      legend.position = "right"
    )
  
  return(p)
}

p_cya <- plot_composition(
  ps_cya,
  "Cyanobacteria baseline composition",
  rank_name = "Genus",
  top_n = 12,
  palette_base = "Paired"
)

p_green <- plot_composition(
  ps_green,
  "Green algae baseline composition",
  rank_name = "Genus",
  top_n = 12,
  palette_base = "Set3"
)

p_euka <- plot_composition(
  ps_euka,
  "Eukaryotic baseline composition",
  rank_name = "Class",
  top_n = 12,
  palette_base = "Spectral"
)

ggsave(
  file.path(fig_dir, "S5_CYA_Baseline_Composition.tiff"),
  p_cya,
  width = 9,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(fig_dir, "S6_GreenAlgae_Baseline_Composition.tiff"),
  p_green,
  width = 9,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(fig_dir, "S7_Euka01_Baseline_Composition.tiff"),
  p_euka,
  width = 9,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# --- 8. Alpha figure ---
p_alpha <- ggplot(alpha_table, aes(x = Marker, y = Observed_ASVs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Observed ASVs",
    title = "Initial biocrust inoculum ASV richness by marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(fig_dir, "S8_AlphaDiversity_T3_Baseline.tiff"),
  p_alpha,
  width = 7,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

# --- 9. Generate manuscript-ready text ---
get_audit_value <- function(marker, column) {
  audit_table %>%
    filter(Marker == marker) %>%
    pull(.data[[column]])
}

cya_asv <- get_audit_value("CYA_16S_Cyanobacteria", "ASVs")
green_asv <- get_audit_value("GreenAlgae_G500_G800", "ASVs")
euka_asv <- get_audit_value("Euka01_18S_Eukaryotes", "ASVs")

cya_reads <- get_audit_value("CYA_16S_Cyanobacteria", "Total_reads")
green_reads <- get_audit_value("GreenAlgae_G500_G800", "Total_reads")
euka_reads <- get_audit_value("Euka01_18S_Eukaryotes", "Total_reads")

baseline_text <- sprintf(
  "### Initial composition of the biocrust inoculum

DNA metabarcoding of the initial T3 inoculum, after LULU curation, contaminant removal, and PCR-replicate filtering, retained %s cyanobacterial 16S ASVs, %s green algal G500/G800 ASVs, and %s broader eukaryotic Euka01 ASVs across the five biological samples from Lachay (S082–S086). Sequencing depth differed among markers, with %s reads for the cyanobacterial marker, %s reads for the green algal marker, and %s reads for Euka01.

The cyanobacterial marker confirmed the presence of a diverse cyanobacterial fraction in the initial inoculum, consistent with the dominance of filamentous cyanobacteria observed under microscopy. The GreenAlgae marker recovered a plastid-based phototrophic signal assigned to chlorophyte and streptophyte lineages, supporting the presence of green algal components in the inoculum. The Euka01 marker recovered a broader eukaryotic assemblage, including phototrophic and heterotrophic lineages. Together, these results confirm that the experimental inoculum represented a mixed biocrust assemblage and provide biological context for interpreting the subsequent functional responses to nutrient addition.

Because metabarcoding was conducted only on the initial inoculum before nutrient addition, these data were used to characterize the starting biological community and were not used to infer treatment-driven changes in microbial community composition.
",
  format(cya_asv, big.mark = ","),
  format(green_asv, big.mark = ","),
  format(euka_asv, big.mark = ","),
  format(cya_reads, big.mark = ","),
  format(green_reads, big.mark = ","),
  format(euka_reads, big.mark = ",")
)

writeLines(
  baseline_text,
  file.path(text_dir, "TEXT_Initial_Biocrust_Inoculum_Results.md")
)

cat("\n======================================================\n")
cat("Baseline metabarcoding summary completed.\n")
cat("Tables saved in:", tab_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
cat("Text saved in:", text_dir, "\n")
cat("======================================================\n")