#!/bin/bash

# ==============================================================================
# Script 00: Automated Extraction, Merging, and ID-Renaming (Eukaryotes 18S)
# Project: Biological Soil Crust Peru (FAIR)
# Purpose: Direct link between Novogene raw data and Clean Sample IDs
# ==============================================================================

# --- 1. SET PATHS  ---
# USER ATTENTION: Replace "/mnt/Your/path/" with the absolute path in your local environment 
# where you downloaded and extracted the root folder "Zenodo_Metabarcoding".
BASE_DIR="/mnt/Your/path/Zenodo_Metabarcoding/01_Data"
RAW_DIR="$BASE_DIR/Raw/Eukaryotes_18S"
METADATA_DIR="$BASE_DIR/Metadata"
METADATA_CSV="$METADATA_DIR/Pierr_metadata.csv"
BARCODE_FILE="$METADATA_DIR/euka01_repbarcodes.txt"
WORK_DIR="$RAW_DIR/seqfiles"

# Create output directory for trimmed files
mkdir -p "$WORK_DIR/euka01trimmed"

# --- 2. DYNAMIC MERGING & AUTOMATIC RENAMING ---
echo "Merging lanes and mapping Novogene names to SampleIDs..."

# Navigate to the raw data folder extracted from Novogene
cd "$WORK_DIR/X204SC23103358-Z01-F001/01.RawData" || { echo "Error: RawData not found"; exit 1; }

for dir in */; do
    # Remove trailing slash for cleaner mapping
    clean_dir=$(echo "$dir" | sed 's/\///g')
    
    # Extract SID using exact match from metadata
    SID=$(grep "^${clean_dir}," "$METADATA_CSV" | cut -d',' -f2 | tr -d '\r')
    
    if [ -n "$SID" ]; then
        echo "Processing: $clean_dir ---> $SID"
        # Merge lanes and move to the working directory
        cat "$clean_dir"/*_1.fq.gz > "$WORK_DIR/${SID}_1.fq.gz"
        cat "$clean_dir"/*_2.fq.gz > "$WORK_DIR/${SID}_2.fq.gz"
    else
        echo "Warning: Folder $clean_dir not found in metadata. Skipping."
    fi
done

cd "$WORK_DIR" || exit
dos2unix "$BARCODE_FILE"

# --- 3. DEMULTIPLEXING & TRIMMING (Main Loop) ---
# Iterating through each merged sample to separate replicates
for r1 in *_1.fq.gz; do
    [ -e "$r1" ] || continue
    SID=${r1/_1.fq.gz}
    echo "--- Processing Sample: $SID ---"

    # Demultiplexing PE reads using sabre
    sabre pe -f "$r1" -r "${SID}_2.fq.gz" -b "$BARCODE_FILE" \
    -u "${SID}_unassigned1.fq" -w "${SID}_unassigned2.fq"

    # Standardize filenames for the trimming step
    mv rep1f.fq "${SID}_rep1f.fq"; mv rep1r.fq "${SID}_rep1r.fq"
    mv rep2f.fq "${SID}_rep2f.fq"; mv rep2r.fq "${SID}_rep2r.fq"
    mv rep3f.fq "${SID}_rep3f.fq"; mv rep3r.fq "${SID}_rep3r.fq"
    mv rep4f.fq "${SID}_rep4f.fq"; mv rep4r.fq "${SID}_rep4r.fq"

    # --- 4. PRIMER TRIMMING (One block per Replicate) ---
    # Optimized using exact primers and reverse complement barcodes
    
    echo "Trimming replicates for $SID..."
    
    # Replicate 1: Rev Barcode GAGTGG
    cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGCCACTC \
             -A GAGTGGCATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA \
             -o "${SID}.rep1.trim1.fq.gz" -p "${SID}.rep1.trim2.fq.gz" \
             "${SID}_rep1f.fq" "${SID}_rep1r.fq" > /dev/null

    # Replicate 2: Rev Barcode CCACGTC
    cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGGACGTGG \
             -A CCACGTCCATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA \
             -o "${SID}.rep2.trim1.fq.gz" -p "${SID}.rep2.trim2.fq.gz" \
             "${SID}_rep2f.fq" "${SID}_rep2r.fq" > /dev/null

    # Replicate 3: Rev Barcode TTCTCAGC
    cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGGCTGAGAA \
             -A TTCTCAGCCATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA \
             -o "${SID}.rep3.trim1.fq.gz" -p "${SID}.rep3.trim2.fq.gz" \
             "${SID}_rep3f.fq" "${SID}_rep3r.fq" > /dev/null

    # Replicate 4: Rev Barcode CCTAG
    cutadapt -a TGGTGCATGGCCGTTCTTAGT...GGTCTGTGATGCCCTTAGATGCTAGG \
             -A CCTAGCATCTAAGGGCATCACAGACC...ACTAAGAACGGCCATGCACCA \
             -o "${SID}.rep4.trim1.fq.gz" -p "${SID}.rep4.trim2.fq.gz" \
             "${SID}_rep4f.fq" "${SID}_rep4r.fq" > /dev/null

    # Cleanup temporary files to save disk space
    rm "${SID}"_rep*.fq "${SID}_1.fq.gz" "${SID}_2.fq.gz"
done

# --- 5. FINAL MOVE ---
# Consolidate all processed files for DADA2 analysis
mv *.trim*.fq.gz ./euka01trimmed/

echo "Step 00 COMPLETE. Data ready for Script 01 (DADA2) in $WORK_DIR/euka01trimmed"