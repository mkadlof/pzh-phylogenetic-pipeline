PZH phylogenetic pipeline
=========================

This project is part of [PleEpiSeq](https://www.pzh.gov.pl/projekty-i-programy/plepiseq/) project.

The pipeline loosely originates from [NextStrain Zika Tutorial](https://github.com/nextstrain/zika-tutorial)

Another project related to PleEpiSeq is [Sequnecing pipline](https://github.com/mkadlof/pzh_pipeline_viral)

To execute bacterial pipeline type

## Minimal Execution on `anat`

Assuming you're in a directory containing the `metadata.txt` file and a `fastas/` directory where each sample’s genome is in a separate FASTA file, run:

```bash
bash nf_pipeline_bacterial_phylo.sh -m metadata.txt \
                                    -i fastas/ \
                                    -t fasta \
                                    -g Salmonella \
                                    -p Salmonella_dummy \
                                    -d /home/michall/git/pzh-phylogenetic-pipeline/ \
                                    --threads 48
```

The `metadata.txt` file must be a tab-separated file with the following **required columns**:

- `strain` – sample identifier, string. 
- `date` – date associated with the sample in `YYYY-MM-DD` format  
- `region` – continent associated with the sample (e.g. `Europe`)  
- `country` – country associated with the sample (e.g. `France`)  
- `division` – sub-country administrative region (e.g. `Ohio`)  
- `city` – city associated with the sample (e.g. `Paris`)  
- `Serovar` – serovar associated with the sample (e.g. `Enteritidis`)  
- `MLST` – MLST ID associated with the sample, number (e.g. `11`)

### Naming Scheme

If the input files in the `fastas/` directory are named like `ERRXYZ.fasta`, the program will convert them to GFF format using Prokka. The resulting GFF files will be named `ERRXYZ_prokka.gff`. As a result, when using FASTA input, the `strain` column in the metadata must match the GFF filenames, i.e., it should be `ERRXYZ_prokka`. If the input type is already GFF (i.e., files are named like `ERRXYZ.gff`), then the `strain` column can simply be `ERRXYZ`.

