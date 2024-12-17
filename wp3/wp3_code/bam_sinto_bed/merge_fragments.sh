#!/bin/bash

# DEFINE PATHS

# Folder to save the bed files in
OUTPUT_BED="/mnt/workspace_stud/stud6/atac/bed"

# Folder for the error files
OUTPUT_ERROR="/mnt/workspace_stud/stud6/atac/error"

# Folder for the merged fragment file
OUTPUT_FRAGMENT="/mnt/workspace_stud/stud6/atac/napkon_fragment"

# Get available threads
NUM_THREADS=$(nproc)

# Path to parallel installtion
parallel_p="/home/stud6/parallel-20241022/src/parallel"


export OUTPUT_BED
export OUTPUT_ERROR
export OUTPUT_FRAGMENT

##################### Concat Files ################################
echo "Concatenating fragment BED files..."

CONCAT_START=$SECONDS
find "${OUTPUT_BED}" -name "*.bed" | $parallel_p -j "$AVAILABLE_PROCS" 'cat {}' > ${OUTPUT_FRAGMENT}/napkon_fragments.bed

echo "Concanted files in $((SECONDS - CONCAT_START)) seconds."
