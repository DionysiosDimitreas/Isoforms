#!/bin/bash

#=============================================================================
# SCRIPT: run_pfam_scan.sh
# DESCRIPTION: Runs hmmscan on a directory of FASTA files.
#
# USAGE: 
#   ./scripts/run_pfam_scan.sh <fasta_directory> <pfam_db_file> <output_directory>
#
# EXAMPLE:
#   ./scripts/run_pfam_scan.sh ./fasta_files ~/db/Pfam-A.hmm ./pfam_results
#=============================================================================

# === Setup ===
# Exit on first error, treat unset variables as an error
set -euo pipefail 

# --- Check for correct number of arguments ---
if [ "$#" -ne 3 ]; then
    echo "ERROR: Incorrect number of arguments."
    echo "Usage: $0 <fasta_directory> <pfam_db_file> <output_directory>"
    exit 1
fi

echo "Job started at $(date)"
echo "Running on host: $(hostname)"

# --- Verify hmmscan exists ---
# (Assumes HMMER is already in the user's system PATH)
which hmmscan > /dev/null || { echo "ERROR: hmmscan not found in PATH. Please install HMMER."; exit 1; }
echo "Found hmmscan:"
hmmscan -h | head -5
echo "-----------------------------------"


# === Paths ===
# Read paths from command-line arguments
FASTA_DIR=$1
PFAM_DB=$2
OUT_DIR=$3
LOG_DIR="${OUT_DIR}/logs"

echo "Input FASTA directory: ${FASTA_DIR}"
echo "Pfam HMM database: ${PFAM_DB}"
echo "Output directory: ${OUT_DIR}"

# --- Create output directories ---
mkdir -p "$OUT_DIR" "$LOG_DIR"


# === Index Pfam database if not already done ===
if [ ! -f "${PFAM_DB}.h3f" ]; then
    echo "Indexing ${PFAM_DB}..."
    hmmpress "${PFAM_DB}"
fi


# === Run hmmscan ===
echo "Starting hmmscan loop..."
for f in "${FASTA_DIR}"/*.fasta; do
    # Get just the filename without the .fasta extension
    base=$(basename "$f" .fasta)
    echo "Processing ${base}..."
    
    hmmscan --cpu 1 --cut_ga \
        --domtblout "${OUT_DIR}/${base}.domtblout" \
        "${PFAM_DB}" "$f" &> "${LOG_DIR}/${base}.log"
done


# === Done ===
echo "-----------------------------------"
echo "Finished all runs."
echo "Results files in: ${OUT_DIR}"
echo "Log files are in: ${LOG_DIR}"
echo "Job finished at $(date)"
