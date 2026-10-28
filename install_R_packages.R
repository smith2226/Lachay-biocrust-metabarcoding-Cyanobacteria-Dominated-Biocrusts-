# Installs the R packages required by the metabarcoding pipeline.
# Versions used for the published analysis are pinned where known; the remaining
# packages are installed at their current release.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))     install.packages("remotes")

# Bioconductor
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "ShortRead"), ask = FALSE)

# CRAN
install.packages(c("dplyr", "tidyr", "stringr", "readr", "magrittr",
                   "ggplot2", "RColorBrewer", "vegan", "here", "Hmisc"))

# LULU is not on CRAN or Bioconductor
remotes::install_github("tobiasgf/lulu")

# Versions used for the published analysis:
#   dada2 1.40.0 | phyloseq 1.56.0 | lulu 0.1.0
cat("Installed. Record your own versions with sessionInfo().\n")
