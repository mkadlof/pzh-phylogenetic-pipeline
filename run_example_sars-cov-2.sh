#!/bin/bash

nextflow run nf_viral_phylogenetic_pipeline.nf \
             --input_fasta data/example_data/sars-cov-2/sars-cov-2.fasta \
             --metadata data/example_data/sars-cov-2/sars-cov-2_metadata.tsv \
             --organism sars-cov-2 \
             -with-dag nf_viral_phylogenetic_pipeline.png \
             -with-trace trace.tsv \
             -with-report report.html \
             -resume

# Remove logs from previous runs
rm \.nextflow.log\.*
