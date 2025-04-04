#!/bin/bash

set -xeuo pipefail

# set -e – przerywa wykonanie skryptu przy pierwszym błędzie.
# set -u – traktuje niezdefiniowane zmienne jako błąd.
# set -o pipefail – powoduje, że błąd w dowolnym elemencie potoku (|) kończy skrypt.
# set -x – włącza tryb debug, czyli wypisuje każdą komendę wraz z podstawionymi wartościami zmiennych.

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
python3 src/filter_low_quality_sequences.py \
              --output_dir "${OUTPUT_DIR}" \
              "${OUTPUT_DIR}/sequence-index.csv"

# Filter sequences
augur filter \
      --sequences "${INPUT_FASTA}" \
      --sequence-index "${OUTPUT_DIR}/sequence-index.csv" \
      --metadata "${METADATA}" \
      --exclude "${OUTPUT_DIR}/invalid_strains.txt" \
      --output "${OUTPUT_DIR}/valid_strains.fasta"

python3 src/find_identical_seqences.py \
                --input "${OUTPUT_DIR}/valid_strains.fasta" \
                --output_dir "${OUTPUT_DIR}"

# Make alignment
augur align \
      --sequences "${OUTPUT_DIR}/valid_strains_unique.fasta" \
      --output "${OUTPUT_DIR}/tree.fasta"

iqtree2 -nt AUTO -s results/tree.fasta -m GTR+G -B 1000 -con -minsup 0.75

# Add missing duplicates
src/insert_missing_dupliated_sequences.py \
      --tree results/tree.fasta.contree \
      --ids results/valid_strains_ident_seq.csv > ${OUTPUT_DIR}/consensus_tree.nwk

#
## Refine tree - currently only for asaignment internal node names
#augur refine \
#      --tree "${OUTPUT_DIR}/tree.nwk" \
#      --output-tree "${OUTPUT_DIR}/refined-tree.nwk" \
#
## Infer ancestral sequences
#augur ancestral \
#    --tree "${OUTPUT_DIR}/refined-tree.nwk" \
#    --alignment "${OUTPUT_DIR}/aligned.fasta" \
#    --output-node-data "${OUTPUT_DIR}/node_data.json" \
