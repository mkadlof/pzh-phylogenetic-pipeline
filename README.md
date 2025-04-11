PZH phylogenetic pipeline
=========================

This project is part of [PleEpiSeq](https://www.pzh.gov.pl/projekty-i-programy/plepiseq/) project.

The pipeline loosely originates from [NextStrain Zika Tutorial](https://github.com/nextstrain/zika-tutorial)

Another project related to PleEpiSeq is [Sequnecing pipline](https://github.com/mkadlof/pzh_pipeline_viral)

To execute bacterial pipeline type

nextflow run nf_pipeline_bacterial_phylo.nf --input_dir /path/to/directory/with/fastas --input_type fasta --results_dir results --threads 20  --main_image pzh_pipeline_phylogenetic --prokka_image staphb/prokka -profile slurm -with-trace -resume
