#!/bin/bash

nextflow run nf_viral_phylogenetic_pipeline.nf \
             --input_fasta data/SARS3.fasta \
             --metadata data/metadata.tsv \
             -with-dag nf_viral_phylogenetic_pipeline.png \
             -with-trace trace.tsv \
             -with-report report.html