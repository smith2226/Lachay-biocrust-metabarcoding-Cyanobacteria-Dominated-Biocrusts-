#!/bin/bash

# ==============================================================================
# Script 00: Green Algae 23S (G500) - Full Processing Pipeline
# Project: Biological Soil Crust Peru (FAIR)
# ==============================================================================

# USER ATTENTION: Replace "/mnt/Your/path/" with the absolute path in your local environment 
# where you downloaded and extracted the root folder "Zenodo_Metabarcoding".
BASE_DIR="/mnt/Your/path/Zenodo_Metabarcoding/01_Data"
RAW_DIR="$BASE_DIR/Raw/Cyano_Green_Mixed"
METADATA_CSV="$BASE_DIR/Metadata/Pierr_metadata.csv"
BARCODE_FILE="$BASE_DIR/Metadata/green_repbarcodes.txt"
WORK_DIR="$RAW_DIR/greenseqfiles"

# Ensure the output directory exists
mkdir -p "$WORK_DIR/Greentrimmed"

# --- 2. AUTO-DETECT SOURCE PATH ---
echo "Detecting RawData location for Algae..."
REAL_RAW_PATH=$(find "$RAW_DIR" -type d -name "01.RawData" | head -n 1)

if [ -z "$REAL_RAW_PATH" ]; then
    echo "Error: 01.RawData not found in $RAW_DIR"
    exit 1
fi

echo "Found RawData at: $REAL_RAW_PATH"

# --- 3. DYNAMIC MERGING ---
echo "Merging lanes for all Green Algae samples..."
cd "$REAL_RAW_PATH" || exit

for dir in */; do
    clean_dir=$(echo "$dir" | sed 's/\///g')
    # Exact SID mapping from metadata
    SID=$(grep "^${clean_dir}," "$METADATA_CSV" | cut -d',' -f2 | tr -d '\r')
    
    if [ -n "$SID" ]; then
        echo "Merging: $clean_dir ---> $SID"
        # Combine all lanes (L1, L2, etc) into merged R1 and R2
        cat "$clean_dir"/*_1.fq.gz > "$WORK_DIR/${SID}_merged_R1.fq.gz"
        cat "$clean_dir"/*_2.fq.gz > "$WORK_DIR/${SID}_merged_R2.fq.gz"
    fi
done

# --- 4. DEMULTIPLEXING & TRIMMING ---
cd "$WORK_DIR" || exit
dos2unix "$BARCODE_FILE"

for r1 in *_merged_R1.fq.gz; do
    [ -e "$r1" ] || continue
    SID=${r1/_merged_R1.fq.gz}
    echo "Processing Sample: $SID"

    # Sabre demux using specific green algae barcodes
    sabre pe -f "$r1" -r "${SID}_merged_R2.fq.gz" -b "$BARCODE_FILE" \
    -u "${SID}_u1.fq" -w "${SID}_u2.fq"

    # Trimming all 4 replicates with G500 Primers
    for r in 1 2 3 4; do
        mv "rep${r}f.fq" "${SID}_rep${r}f.fq"
        mv "rep${r}r.fq" "${SID}_rep${r}r.fq"
        
        echo "Trimming $SID Rep $r..."
        # Using the specific G500 23S primers
        cutadapt -a GAATGAGTACAATCTAAACCCCTTAAC...CCTGTTGGTCTGTAGGACCGGAGTAATG \
                 -A CATTACTCCGGTCCTACAGACCAACAGG...GTTAAGGGGTTTAGATTGTACTCATTC \
                 -o "${SID}.rep${r}.trim1.fq.gz" -p "${SID}.rep${r}.trim2.fq.gz" \
                 "${SID}_rep${r}f.fq" "${SID}_rep${r}r.fq" > /dev/null
    done
    
    # Clean up merged and temporary rep files to save space
    rm "${SID}_merged"* "${SID}_rep"* "${SID}_u"*
done

# --- 5. FINAL MOVE ---
mv *.trim*.fq.gz ./Greentrimmed/
echo "Step 00 COMPLETE. Green Algae files ready in $WORK_DIR/Greentrimmed"