PZH phylogenetic pipeline
=========================

This project is part of [PleEpiSeq](https://www.pzh.gov.pl/projekty-i-programy/plepiseq/) project.

The pipeline loosely originates from [NextStrain Zika Tutorial](https://github.com/nextstrain/zika-tutorial)

Another project related to PleEpiSeq is [Sequnecing pipline](https://github.com/mkadlof/pzh_pipeline_viral)

To execute bacterial pipeline type

## Minimal Execution on `anat`

Assuming you're in a directory containing the `metadata_new.txt` file and a `fastas/` directory where each sample’s genome is in a separate FASTA file, run:

```bash
bash nf_pipeline_bacterial_phylo.sh -m metadata_new.txt -i fastas/ -t fasta -g Salmonella -p Salmonella_dummy -d /home/michall/git/pzh-phylogenetic-pipeline/ --threads 48
```
