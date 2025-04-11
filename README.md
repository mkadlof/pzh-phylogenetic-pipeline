PZH phylogenetic pipeline
=========================

This project is part of [PleEpiSeq](https://www.pzh.gov.pl/projekty-i-programy/plepiseq/) project.

The pipeline loosely originates from [NextStrain Zika Tutorial](https://github.com/nextstrain/zika-tutorial)

Another project related to PleEpiSeq is [Sequnecing pipline](https://github.com/mkadlof/pzh_pipeline_viral)

To execute bacterial pipeline type

nextflow run nf_pipeline_bacterial_phylo.nf --input_dir /path/to/dir/with/fastas --input_type fasta --results_dir ./results --prokka_image staphb/prokka --main_image pzh_pipeline_phylogenetic --metadata /mnt/raid/michall/test_bakterie_filogeneryka_2_pipeline/metadata.txt --threads 30 -with-trace -resume -profile slurm
