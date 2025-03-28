#!/bin/bash

INPUT_FASTA=$1
METADATA=$2
OUTPUT_DIR=$3


# if output directory does not exist, create it
if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir -p ${OUTPUT_DIR}
fi

augur index \
      --sequences ${INPUT_FASTA} \
      --output ${OUTPUT_DIR}/sequence-index.csv

python3 src/filter_low_quality_sequences.py ${OUTPUT_DIR}/sequence-index.csv --output_dir ${OUTPUT_DIR}

augur filter \
      --sequences "${INPUT_FASTA}" \
      --sequence-index "results/sequence-index.csv" \
      --metadata "${METADATA}" \
      --group-by month \
      --sequences-per-group 3 \
      --exclude "${OUTPUT_DIR}/invalid_strains.txt" \
      --output "${OUTPUT_DIR}/valid_strains.fasta"

augur align \
      --sequences "${OUTPUT_DIR}/valid_strains.fasta" \
      --output "${OUTPUT_DIR}/aligned.fasta" \

augur tree \
      --alignment "${OUTPUT_DIR}/aligned.fasta" \
      --output "${OUTPUT_DIR}/tree.nwk" \
      --method "iqtree"
