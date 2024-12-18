#!/bin/bash

# DEFINE PATHS

# Path to the file with the bam file paths
INPUT_BAM_FILES="/mnt/workspace_stud/stud6/pul52_napkon_atac_bams.txt"
#test_data.txt"
#pul52_napkon_atac_bams.txt"

# Folder to save the bed files in
OUTPUT_BED="/mnt/workspace_stud/stud6/atac/bed"

# Folder to store the modified BAM files (adding CELL BARCODE)
OUTPUT_BAM_BC="/mnt/workspace_stud/stud6/atac/bam_cb"

# Folder for the error files
OUTPUT_ERROR="/mnt/workspace_stud/stud6/atac/error"

# Folder for the merged fragment file
OUTPUT_FRAGMENT="/mnt/workspace_stud/stud6/atac/napkon_fragment"

# Get available threads
NUM_THREADS=$(nproc)

# Path to parallel installtion
parallel_p="/mnt/workspace_stud/stud6/parallel-20241022/src/parallel"

# Create output dir if not exists
mkdir -p "$OUTPUT_BED"
mkdir -p "$OUTPUT_BAM_BC"
mkdir -p "$OUTPUT_ERROR"
mkdir -p "$OUTPUT_FRAGMENT"

# Create error logging file
ERROR_LOG="$OUTPUT_ERROR/error2.log"
> "$ERROR_LOG"

# Check if sinto, samtools and gnu parallel are installed and avail
if ! command -v sinto &> /dev/null
then
    echo "sinto installation not found."
    exit 1
fi

if ! command -v $parallel_p &> /dev/null
then
    echo "GNU Parallel installation not found. Parallel processing will not be available." >> "$ERROR_LOG"
    exit 1
fi

if ! command -v samtools &> /dev/null
then
    echo "samtools installation not found."
    exit 1
fi


##################### Define Functions ###########################
add_barcodes() {
    # Get arg which should contain bam file
    BAM_FILE=$1
    # Create output filenames
    BASE_NAME=$(basename "$BAM_FILE" .bam)
    
    # Extract the patient ID or cell barcode
    #CELL_BARCODE=$(echo "$BASE_NAME" | awk -F'[_.]' '{for(i=1;i<=NF;i++) if($i ~ /^Plate/) {print $i"_"$(i+1); exit}}')
    CELL_BARCODE="$BASE_NAME"
    echo "CELL_BARCODE: $CELL_BARCODE"
    echo "Processing $BAM_FILE ..."
    START_TIME=$SECONDS
    
    # Check if BAM_FILE exists
    if [ ! -s "$BAM_FILE" ]; then
        echo "Error: $BAM_FILE is empty or missing." >> "$ERROR_LOG"
        return 1
    fi      
    
    # Add cell barcodesS
    OUTFILE="${OUTPUT_BAM_BC}/${BASE_NAME}.cb.bam"
    
    # Check if BAM_FILE exists
    if [ ! -s "$OUTFILE" ]; then
    
      echo "Adding barcode to $BAM_FILE..."
    
      # Create index file in OUTPUT_BAM_BC folder
      if [ ! -f "${OUTPUT_BAM_BC}/${BASE_NAME}.bai" ]; then
        echo "Indexing $BAM_FILE..."
        STEP_START=$SECONDS
        if ! samtools index --threads "$MAX_PROCS_PER_FILE" "$BAM_FILE" "${OUTPUT_BAM_BC}/${BASE_NAME}.bai"; then
          echo "Error: Failed to index $BAM_FILE." >> "$ERROR_LOG"
          return 1
        fi
      fi
      
      echo "Indexing completed in $((SECONDS - STEP_START_INDEX)) seconds."
      STEP_START=$SECONDS
    
      if ! python "/home/stud6/wp3/wp3_code/bam_sinto_bed/add_barcode.py" "$BAM_FILE" "$CELL_BARCODE" "$OUTFILE"; then
          echo "Error: Failed to add barcode to $BAM_FILE." >> "$ERROR_LOG"
      fi
      echo "Barcode added in $((SECONDS - STEP_START)) seconds."
  
      echo "Removing index file: ${OUTPUT_BAM_BC}/${BASE_NAME}.bai"
      
      # Remove the index file because we need a new one to add fragments with the barcode column
      rm -f "${OUTPUT_BAM_BC}/${BASE_NAME}.bai"
    fi
    
    
    #echo "MV bc.bam file to $BAM_FILE"
    #mv ${BAM_FILE%.bam}.bc.bam $BAM_FILE
}

export -f add_barcodes

# Function to process the bam files
fragment_bam() {
    # Get arg which should contain bam file
    BAM_FILE=$1
    # Create output filenames
    BASE_NAME=$(basename "$BAM_FILE" .cb.bam)
    
    # STUDY: Is compression with gz faster to process?
    FRAGMENTS_FILE="$OUTPUT_BED/$BASE_NAME.bed"
    
    #TN5_ADJUSTMENT="+4/-5" # Default is fine according to Jan, thus, do not parse arg
    MAX_TN5_DISTANCE=1000   # Default of PeakQC
    MIN_MAPPING_QUALITY=30 # According to Jan
    
    if [ ! -s "$FRAGMENTS_FILE" ]; then
    
      # Files are sorted by coordinate 
      if [ ! -f "${BAM_FILE}.bai" ]; then
        echo "Indexing $BAM_FILE..."
        STEP_START=$SECONDS
        if ! samtools index --threads "$MAX_PROCS_PER_FILE" "$BAM_FILE"; then
          echo "Error: Failed to index $BAM_FILE." >> "$ERROR_LOG"
          return 1
        fi
      fi
          
      
      # Fragment file
      echo "Creating fragment file..."
      STEP_START=$SECONDS
      if ! sinto fragments -b "$BAM_FILE" -f "$FRAGMENTS_FILE" -p "$MAX_PROCS_PER_FILE" --max_distance "$MAX_TN5_DISTANCE" --min_mapq "$MIN_MAPPING_QUALITY"; then
          echo "Error: Failed to create fragment file for $BAM_FILE." >> "$ERROR_LOG"
          return 1
      fi
          
      # Remove the index file to save space on disk
      rm -f "$BAM_FILE.bai"
    fi
    
    echo "Processed $BAM_FILE in $((SECONDS - START_TIME)) seconds."
}
# Make avail to GNU parallel
export -f fragment_bam
export OUTPUT_BED
export OUTPUT_BAM_BC
export OUTPUT_ERROR
export OUTPUT_FRAGMENT

###################################################################
start_time=$(date +%s)

# sinto needs sorted and indexed bam files. We do not have index files but the 
# files are sorted by coordinate

# Hwo many processes, I/o bandwith and memory avail?
AVAILABLE_PROCS=$(nproc)
# the max processes per file are 24 because we have 24 chromosome and one process will handle one chromosome
# check the scaling behavior of sinto but also index
MAX_PROCS_PER_FILE=22


# Check how many processes are avail 
if [ "$AVAILABLE_PROCS" -le "$MAX_PROCS_PER_FILE" ]; then
    MAX_PROCS_PER_FILE=$AVAILABLE_PROCS
fi

export MAX_PROCS_PER_FILE

echo "Using $MAX_PROCS_PER_FILE processors for each BAM file."

#BAM_DIR="/mnt/workspace_stud/stud6/atac/bam/"
#BAM_FILES=("$BAM_DIR"/*.bam)

mapfile -t BAM_FILES < "$INPUT_BAM_FILES"


FILES_COUNT=${#BAM_FILES[@]}


SPLIT=$(( FILES_COUNT / 2 ))

SUBSET_BAM=("${BAM_FILES[@]:SPLIT:FILES_COUNT}")

echo "Number of BAM files in subset to process: ${#SUBSET_BAM[@]}"

##################### ADDING BARCODES #############################
echo "Adding the cell barcode tags..."

$parallel_p -j 10 add_barcodes ::: "${SUBSET_BAM[@]}"

##################### Fragment BAM ################################
# Process multiple bam files if nprocs > max_procs

BAM_CB_FILES=("$OUTPUT_BAM_BC"/*.cb.bam)

SUBSET_BAM_CB=("${BAM_CB_FILES[@]:SPLIT:FILES_COUNT}")

FILES_IN_PARALLEL=$((AVAILABLE_PROCS / MAX_PROCS_PER_FILE))

echo "Fragment ${#SUBSET_BAM_CB[@]} files..."

if [ "$FILES_IN_PARALLEL" -gt "${#SUBSET_BAM_CB[@]}" ]; then
    FILES_IN_PARALLEL=$FILES_COUNT
fi


if [ "$AVAILABLE_PROCS" -gt "$MAX_PROCS_PER_FILE" ]; then
    echo "More than $MAX_PROCS_PER_FILE processors available. Processing ${FILES_IN_PARALLEL} BAM files in parallel."

    # Process the BAM files in parallel, each using MAX_PROCS_PER_FILE processes
    $parallel_p -j $FILES_IN_PARALLEL fragment_bam ::: "${SUBSET_BAM_CB[@]}"

else
    echo "$AVAILABLE_PROCS processors available. Processing one BAM file at a time."

    # Process each BAM file one by one using the available processors (>=22)
    for BAM_FILE in "${SUBSET_BAM_CB[@]}"; do
        fragment_bam "$BAM_FILE";
    done
fi

###################################################################
echo "All BAM files have been processed."

end_time=$(date +%s)
run_time=$((end_time - start_time))
echo "Total run time: $run_time seconds."
