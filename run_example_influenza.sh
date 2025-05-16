#!/bin/bash

nextflow run nf_viral_phylogenetic_pipeline.nf \
             --input_fasta data/example_data/influenza/ \
             --metadata data/example_data/influenza/influenza_metadata.tsv \
             --organism influenza \
             --threshold_Ns 0.1 \
             --threshold_ambiguities 1.0 \
             -with-dag nf_viral_phylogenetic_pipeline.png \
             -with-trace trace.tsv \
             -with-report report.html \
             -resume

# Remove logs from previous runs
rm \.nextflow.log\.*
