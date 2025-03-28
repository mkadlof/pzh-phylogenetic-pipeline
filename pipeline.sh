#!/bin/bash

usage() {
    echo "Usage: $0 -i INPUT_FASTA -m METADATA -o OUTPUT_DIR"
    echo "  -i INPUT_FASTA   Path to the input FASTA file"
    echo "  -m METADATA      Path to the metadata file"
    echo "  -o OUTPUT_DIR    Path to the output directory"
    exit 1
}

while getopts ":i:m:o:" opt; do
    case ${opt} in
        i )
            INPUT_FASTA=$OPTARG
            ;;
        m )
            METADATA=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        \? )
            usage
            ;;
    esac
done

# Check if all required arguments are provided
if [ -z "${INPUT_FASTA}" ] || [ -z "${METADATA}" ] || [ -z "${OUTPUT_DIR}" ]; then
    usage
fi

# if output directory does not exist, create it
if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi

# Index sequences
augur index \
      --sequences "${INPUT_FASTA}" \
      --output "${OUTPUT_DIR}/sequence-index.csv"

# Identify low quality sequences
python3 src/filter_low_quality_sequences.py "${OUTPUT_DIR}/sequence-index.csv" --output_dir "${OUTPUT_DIR}"

# Filter sequences
augur filter \
      --sequences "${INPUT_FASTA}" \
      --sequence-index "${OUTPUT_DIR}/sequence-index.csv" \
      --metadata "${METADATA}" \
      --group-by month \
      --sequences-per-group 3 \
      --exclude "${OUTPUT_DIR}/invalid_strains.txt" \
      --output "${OUTPUT_DIR}/valid_strains.fasta"

# Make alignment
augur align \
      --sequences "${OUTPUT_DIR}/valid_strains.fasta" \
      --output "${OUTPUT_DIR}/aligned.fasta"

# Build tree
augur tree \
      --alignment "${OUTPUT_DIR}/aligned.fasta" \
      --output "${OUTPUT_DIR}/tree.nwk" \
      --method "iqtree"