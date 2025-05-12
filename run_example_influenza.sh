#!/bin/bash

nextflow run nf_viral_phylogenetic_pipeline.nf \
             --input_fasta data/example_data/influenza/ \
             --metadata data/example_data/influenza/influenza_metadata.tsv \
             --organism influenza \
             -with-dag nf_viral_phylogenetic_pipeline.png \
             -with-trace trace.tsv \
             -with-report report.html \
             -resume

# Remove logs from previous runs
rm \.nextflow.log\.*
