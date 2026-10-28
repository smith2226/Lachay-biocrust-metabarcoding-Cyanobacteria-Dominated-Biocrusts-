# Metabarcoding pipeline — Lachay biocrust (Soil Systems)

Bioinformatic pipeline used to process the DNA metabarcoding data reported in:

> Cisneros, P.; Dunshea, G.; Martin, M.D.; Arana, C. **Beyond Phosphorus: Nutrient Medium
> Composition Shapes Early Growth in Cyanobacteria-Dominated Biocrusts from a Peruvian Fog
> Oasis.** *Soil Systems* **2026**. https://doi.org/__xxxxxx__

This repository contains **only the metabarcoding pipeline**: demultiplexing and primer
trimming, DADA2 sequence inference, LULU curation, contaminant removal, phyloseq assembly
and taxonomic assignment for the three markers, plus the script that summarises the initial
inoculum community reported in Section 3.1.

Everything else associated with the article is distributed elsewhere and is **not duplicated
here**:

| Material | Where it lives |
|---|---|
| Microcosm mass-gain, chlorophyll *a*, net CO₂ exchange and nutrient-input analyses | Supplementary Code S1 (published with the article) |
| Curated ASV tables, taxonomic assignments, sample metadata and controls | Supplementary Data S4 and S5 |
| Curated phyloseq objects for this study | Zenodo [10.5281/zenodo.21655411](https://doi.org/10.5281/zenodo.21655411) |
| Raw reads and full pipeline outputs, 16S + 23S | Zenodo [10.5281/zenodo.21654921](https://doi.org/10.5281/zenodo.21654921) |
| Raw reads and full pipeline outputs, 18S | Zenodo [10.5281/zenodo.21655235](https://doi.org/10.5281/zenodo.21655235) |
| Laboratory protocol for the metabarcoding workflow | Supplementary Methods S1 |

## Pipeline

Each marker follows the same sequence of steps. Scripts are numbered in execution order.

| Step | Tool | What it does |
|---|---|---|
| `00_*_Demux_Trim.sh` | sabre + cutadapt | Lane merging, demultiplexing by PCR replicate, primer trimming |
| `01_*_DADA2_Inference.R` | DADA2 | ASV inference with a NovaSeq-corrected error model |
| `02_*_LULU_Curation.R` | LULU + BLAST+ | Collapses redundant ASVs and removes polymerase artefacts |
| `03_*_Contamination_Filtering.R` | phyloseq | Removes contaminants identified in field blanks, extraction blanks and PCR negative controls; enforces the rule that an ASV must be detected in at least two of four PCR replicates |
| `04_*_Phyloseq_Export.R` | phyloseq | Builds the master phyloseq objects with sample metadata |
| `05_*_TaxonomicAssignment*.R` | DADA2 `assignTaxonomy` (`minBoot = 80`) | CyanoSeq v1.3 for the cyanobacterial 16S marker; PR2 v5.1.0 for the plastid 23S and eukaryotic 18S markers |
| `06_*_Subsetting.R` | phyloseq | Splits the master objects into taxonomic groups |
| `07_*_DistributingDatasetsFolders.R` | phyloseq | Writes the per-study subsets |

```
02_Analysis_Core/
├── scripts/
│   ├── 01_Eukaryotes_18S/       Euka01 marker
│   ├── 02_Cyanobacteria_16S/    Cyanobacterial marker
│   └── 03_GreenAlgae_23S/       Plastid G500/G800 marker
├── objects/                     written by the pipeline
└── results/                     written by the pipeline
03_Paper01_Lomas_Biomass/
└── scripts/metabarcoding/01_Baseline_Community.R
```

The directory numbering mirrors the project's working structure, which the scripts'
relative paths assume. Empty `objects/` and `results/` folders are included so the pipeline
can write into them on a fresh clone.

`03_Paper01_Lomas_Biomass/scripts/metabarcoding/01_Baseline_Community.R` reads the curated phyloseq objects and
produces the ASV counts, sequencing depths and taxonomic summaries reported in Section 3.1
and Supplementary Table S3. It also writes exploratory composition and alpha-diversity
figures that were not included in the article.

## Software

| | Version |
|---|---|
| R | 4.5.1 |
| DADA2 | 1.40.0 |
| LULU | 0.1.0 ([tobiasgf/lulu](https://github.com/tobiasgf/lulu)) |
| phyloseq | 1.56.0 |
| sabre | 1.000 |
| cutadapt | 4.7|
| BLAST+ | 2.16.0+ |

Reference databases: **CyanoSeq v1.3** (cyanobacterial 16S) and **PR2 v5.1.0** (plastid 23S
and eukaryotic 18S).

Install the R stack with:

```bash
Rscript install_R_packages.R
```

The command-line tools (`sabre`, `cutadapt`, `blastn`, `makeblastdb`) are listed in
`environment.yml` and can be installed with `conda env create -f environment.yml`.

## Running the pipeline

The scripts expect the raw data and reference databases alongside this repository. Download
the two raw-data records from Zenodo and reconstruct this layout:

```bash
git clone <repository URL> GitHub_Metabarcoding
```

The local folder **must be named `GitHub_Metabarcoding`**: the `here()` calls inside the R
scripts resolve against it. Then download the two raw-data records from Zenodo and build
this layout around it:

```
project_root/
├── GitHub_Metabarcoding/         this repository
└── Zenodo_Metabarcoding/
    └── 01_Data/
        ├── Raw/
        │   ├── Cyano_Green_Mixed/
        │   │   ├── cyanoseqfiles/Cyanotrimmed/      from 10.5281/zenodo.21654921
        │   │   └── greenseqfiles/Greentrimmed/     from 10.5281/zenodo.21654921
        │   └── Eukaryotes_18S/seqfiles/euka01trimmed/   from 10.5281/zenodo.21655235
        ├── Metadata/Pierr_metadata.csv
        └── Databases/
            ├── PR2/pr2_version_5.1.0_SSU_dada2.fasta.gz
            └── Cyanoseq/CyanoSeqV1.3_dada2.fastq.gz
```

Then set the data directory at the top of each `00_*.sh` script (the placeholder reads
`/mnt/Your/path/`) and run the numbered scripts in order for each marker, with
`project_root/` as the R working directory. Approximately 80 GB of free disk space is needed for the
raw data. A Linux environment (or WSL2) is assumed for the shell scripts.

To reproduce only the Section 3.1 summaries without downloading raw reads (~20 MB), take the
curated phyloseq objects from [10.5281/zenodo.21655411](https://doi.org/10.5281/zenodo.21655411)
or from Supplementary Data S4 and run `scripts/04_Baseline_Community/01_Baseline_Community.R`.

## Scope

Metabarcoding data in this study characterise the **initial biocrust inoculum only**, before
nutrient treatment. They are not used to infer treatment-driven changes in community
composition.

Raw demultiplexed reads for the article-specific samples and their controls form part of an
ongoing larger study and are not yet publicly available. Requests should be directed to the
corresponding author (caranab@unmsm.edu.pe).

## License

Code is released under the MIT License (see `LICENSE`).
