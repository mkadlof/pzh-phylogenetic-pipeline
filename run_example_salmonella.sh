#!/bin/bash

# Example: Run phylogentic pipline on Salmonella samples
# All genomes are stored in a single directory in distinct FASTA files.
# Nextflow will use LOCAL to manage process execution.
# External databases are NOT part of this repo. The default below is the shared
# NFS resource; override it if your databases live somewhere else.

PATH_TO_EXTERNAL_DATABASES="/mnt/unity_nfs/external_databases"

bash nf_pipeline_bacterial_phylo.sh --metadata data/example_data/salmonella/metadata_salmonella.txt \
                                    --inputDir data/example_data/salmonella/fastas \
				    --genus Salmonella \
				    --inputType fasta \
				    --projectDir `pwd` \
				    -x "local" \
				    --db "${PATH_TO_EXTERNAL_DATABASES}" \
				    --results_prefix test_salmonella \
				    --results_dir results