PZH phylogenetic pipeline
=========================

This project is part of [PleEpiSeq](https://www.pzh.gov.pl/projekty-i-programy/plepiseq/) project, co-funded by the European Union.

This repository contains a Nextflow pipeline for phylogenetic analysis of viral and bacterial genomes. The pipeline is designed to be modular and can be easily extended to include additional steps or tools.

Pipeline overview
-----------------
Following links provide to PNG files with DAG overview of the pipeline:
- [Sars pipeline](dag_png/nf_sars_phylogenetic_pipeline.png)
- [Influenza pipeline](dag_png/nf_influenza_phylogenetic_pipeline.png)
- [Bacterial pipeline](dag_png/nf_bacterial_phylogenetic_pipeline.png)

Quick start
-----------

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html)
2. Install [Docker](https://docs.docker.com/engine/install/)
3. Clone this repository
4. Build docker pipeline images:
   ```bash
   docker build -t pzh_pipeline-phylo -f Dockerfiles/Dockerfile .
   ```
5. Pull [staphb/prokka image](https://hub.docker.com/r/staphb/prokka) image from DockerHub using:
   ```bash
   docker pull staphb/prokka:latest
   ```
6. Set up your Nextflow configuration file. Recommended way is to copy the example configuration file in the repository root.
   ```bash
   cp nextflow.config.template nextflow.config
   ```
   **Important**: adjust the content of this file to your particular environment and needs.
7. [Optional] Run example scripts to test the pipeline:
   ```bash
   ./run_example_sars-cov-2.sh
   ./run_example_influenza.sh
   ./run_example_salmonella.sh
   ./run_example_rsv.sh   
    ```
   Expected outputs of these scripts can also be found under: `data/example_data/end-2-end/`. 

   Note: `run_example_influenza.sh` and `run_example_rsv.sh` use `slurm` executor, modify these files if you want to use `local` executor. `run_example_salmonella.sh` points `PATH_TO_EXTERNAL_DATABASES` at the shared `/mnt/unity_nfs/external_databases`; adjust it if your databases live elsewhere.
   

8. [Optional] Validate example outputs with pytest. Default outputs of the example scripts are in `results/test_*` directories.
   ```
   pytest tests/end-2-end/test_end-2-end.py --data-dir result/test_sars
   pytest tests/end-2-end/test_end-2-end.py --data-dir result/test_influenza
   pytest tests/end-2-end/test_end-2-end.py --data-dir result/test_rsv
   pytest tests/end-2-end/test_end-2-end.py --data-dir result/test_salmonella
   ```
   or using pre-prepared goldens (e.g. for sars-cov-2):

   ```
   pytest tests/end-2-end/test_end-2-end.py --data-dir data/example_data/end-2-end/result/test_sars
   ```

   Dependencies: pytest, ete3, fastjsonschema, requests

Related projects
----------------

The pipeline loosely originates from [NextStrain Zika Tutorial](https://github.com/nextstrain/zika-tutorial)

Another project related to PleEpiSeq is [Sequnecing pipline](https://github.com/mkadlof/pzh_pipeline_viral)

-----------------------------------------------------------------

# Viral Phylogenetic Pipeline

This pipeline constructs viral phylogenies using raw FASTA files only.

## Input

The phylogenetic pipeline requires two types of input:
1. A tab-separated metadata file (see example below).
2. A directory containing full genomic sequences of analyzed samples in FASTA format.  

### Input directory example
Each genome must be stored in a **separate gzipped FASTA file (`.fasta.gz`)**.  
For multi-segment viruses (e.g. *Influenza*), all genomic segments must be included in a **single FASTA file per sample**.

Example directory layout:
```
/some/path/
├── Sample1.fasta.gz
├── Sample2.fasta.gz
└── Sample3.fasta.gz
```

See `data/example_data/` for reference examples.

### FASTA headers
FASTA headers **must** follow the standardized format established in the WGS pipeline:
> Segment_name|Sample_name

Where: 
- `segment_name` – identifier of the genomic segment (e.g., `chr1_PB2`, `MN908947.3`).  
  For single-segment viruses such as *SARS-CoV-2* or *RSV*, this usually corresponds to the reference sequence identifier.
- `sample_name` – unique sample identifier; must also appear in the metadata file under the `strain` column.

For example:
```
>chr4_HA|Influenza_sample2
ATGGAGAAA...
>chr6_NA|Influenza_sample2
ATAAAAGCC
```

### Metadata linkage
The `sample_name` must appear in the metadata file (in the `strain` column).  
Only this portion (after the `|` character) is matched against the metadata.


## Minimal Execution

1. Create a working directory where you want to store the results.
2. Copy the `nf_pipeline_viral_phylo.sh` script from the repository’s root directory into your working directory.
3. (Optional) Copy valid metadata and fasta file from `/data/example_data/SPECIES/` into the working directory. Use data from viral species.

Call the wrapper script

```bash
bash nf_pipeline_viral_phylo.sh --inputDir PATH_TO_DIRECTORY_WITH_FASTAS \ 
                                --metadata PATH_TO_METADATA_FILE 
                                --organism SELECTED_SPECIES 
                                --results_prefix PROJECT_NAME 
                                --projectDir PATH_TO_REPOSITORY
```

e.g. if one copy data from `/data/example_data/sars-cov-2` to a working directory, and cloned this repo to `/home/my_user/plepiseq-phylogenetic-pipeline`, the command would be

```bash
bash nf_pipeline_viral_phylo.sh --inputDir sars-cov-2 \
                                --metadata sars-cov-2/sars-cov-2_metadata.tsv \
                                --organism sars-cov-2 \
                                --results_prefix sars_example \
                                --projectDir /home/my_user/plepiseq-phylogenetic-pipeline
```

To see all available options and customize your run, use:

```bash
bash nf_pipeline_viral_phylo.sh -h
```

## Metadata Format

The metadata file must be a tab-separated file with the following required columns:

- `strain` – Sample identifier (must match the filename, excluding extension). Sample identifier must be a part of a header in FASTA file
- `virus` - Virus name (sars-cov-2, influenza, or rsv)
- `date` – Collection date in `YYYY-MM-DD` format
- `country` – Country name (e.g., `France`)
- `city` – City name (e.g., `Paris`)
- `type` – Additional classification column representing suptype for Influenza (e.g `H1N1` )or type for for RSV (e.g. `A`). For SARS-CoV-2 can be identical with `virus` column

Other columns are optional and can be used for additional metadata.

To ensure homogeneity of input data, the pipeline applies **safeguards**.  
The default safeguard level is **`type`**, meaning that all samples in one run must share the same `type` value.

Depending on the configured safeguard level, the pipeline will not execute if:

- More than one unique **`virus`** is present in the `virus` column.
- More than one unique **`type`** is present in the `type` column (**default safeguard**).
- The column selected by `--map_detail` (`city` by default, or `country`) contains empty values. Empty join keys duplicate rows when metadata is merged with coordinates.
- The `date` column contains a value that is not `YYYY-MM-DD`, or every sample shares the exact same date and `--clockrate` was not provided. TimeTree cannot estimate a clock rate without variation in sampling dates; either add samples from another date or pass `--clockrate` explicitly.

- This prevents accidental mixing of heterogeneous datasets (e.g., different viruses, types, or geographic origins) in a single phylogenetic analysis run.
---

-----------------------------------------------------------------

# Bacterial Phylogenetic Pipeline

This pipeline constructs bacterial phylogenies using annotated assemblies or raw FASTA files.

## Input

Like the viral counterpart the bacterial phylogenetic pipeline also requires two types of input:
1. A tab-separated metadata file (see example below).
2. A directory containing genomic assemblies of analyzed isolates in FASTA format.
3. A path to external databases. Consult MST Tree section below.
### Input directory example
Each isolate must be represented by a **single gzipped FASTA file (`.fasta.gz`)** containing all assembled contigs for that sample.
Example directory layout:

```
/some/path/
├── Sample1.fasta.gz
├── Sample2.fasta.gz
└── Sample3.fasta.gz
```

### FASTA headers
Unlike viral genomes, the **FASTA headers are not relevant** for bacterial input.  
Only the **file names** are used to identify samples and link them to metadata.  
Each `.fasta.gz` file should therefore include all contigs for a single isolate.

> **Tip:** The filenames (without the `.fasta.gz` extension) must match the `strain` column in the metadata file.

## Minimal Execution

Follow these steps to run the pipeline with minimal setup:

1. Create a working directory where you want to store the results.
2. Copy the `nf_pipeline_bacterial_phylo.sh` script from the repository’s root directory into your working directory.
3. (Optional) Copy a valid metadata file e.g. `metadata_salmonella.txt` file and a `fastas/` directory with uncompressed genome FASTA files into the working directory.  Example files are available in the `data/example_data/salmonella` directory of the repository.

Assuming your working directory contains `metadata_salmonella.txt` and a `fastas/` directory, cloned this repo to `/home/my_user/plepiseq-phylogenetic-pipeline`, and you’ve built/pulled the required Docker images as described in Quick Start section, and external databases are located in `/mnt/unity_nfs/external_databases` run:

```bash
bash nf_pipeline_bacterial_phylo.sh --metadata metadata_salmonella.txt \
                                    --inputDir fastas/ \
                                    --inputType fasta \
                                    --genus Salmonella \
                                    --projectDir /home/my_user/plepiseq-phylogenetic-pipeline \
                                    --results_prefix Salmonella_test \
                                    --db /mnt/unity_nfs/external_databases
```

`--db` defaults to `/mnt/unity_nfs/external_databases`, the shared NFS resource mounted on all compute nodes, so it can be omitted in that environment. Pass it explicitly if your databases live elsewhere.

To see all available options and customize your run, use:

```bash
bash nf_pipeline_bacterial_phylo.sh -h
```

One can also provide paths to metadata and directory with fasta files. There is no need to copy them to working directory.

---

## Metadata Format

The metadata file must be a tab-separated file with the following required columns:

- `strain` – Sample identifier (must match the filename, excluding extension)
- `date` – Collection date in `YYYY-MM-DD` format
- `region` – Continent (e.g., `Europe`)
- `country` – Country name (e.g., `France`)
- `division` – Sub-national region (e.g., `Ohio`)
- `city` – City name (e.g., `Paris`)
- `Serovar` – Serovar designation (e.g., `Enteritidis`)
- `MLST` – MLST ID (e.g., `11`)
- `cgMLST` – cgMLST ID (e.g., `110234`)
- `HC5` – HierCC cluster ID at threshold 5 (e.g., `13`)
- `HC10` – HierCC cluster ID at threshold 10 (e.g., `13`)

The pipeline includes strict safeguards to ensure homogeneity of input data. The pipeline will not execute if:

- Different **serotypes** (`Serovar` column in the metadata file) are provided together.
- The column selected by `--map_detail` (`city` by default, or `country`) is missing or contains empty values. Empty join keys duplicate rows when metadata is merged with coordinates.
- The `date` column contains a value that is not `YYYY-MM-DD`, or every sample shares the exact same date and `--clockrate` was not provided. TimeTree cannot estimate a clock rate without variation in sampling dates; either add samples from another date or pass `--clockrate` explicitly.

---

## Input File Naming

Input files must be provided in either **FASTA** or **GFF** format.

- If FASTA files are used, **Prokka** will be automatically invoked to annotate them to GFF format.
- The `strain` value in the metadata must **exactly match** the filename (excluding extension).

**Example**:  
For a file named `ERRXYZ.fasta`, the corresponding `strain` value in the metadata file must be `ERRXYZ`.


-----------

# WGS2phylo

A helper script to prepare metadata file based on the results of our [Sequnecing pipline](https://github.com/mkadlof/plepiseq-wgs-pipeline). 
With `--with-fasta` a phylogenetic pipeline-ready fasta input can also be prepared. The generated FASTA files and metadata are directly compatible with the viral and bacterial phylogenetic pipelines described above.
Use `--without-fasta` to skip fasta file processing.

## Example data 
- [WGS output for bacterial](data/example_data/WGS2phylo/)

## Execution

### Campylobacter 
```
python3 bin/WGS2Phylo.py --output_dir data/example_data/WGS2phylo/results_Campylobacter/ --organism campylobacter --supplemental-file data/example_data/WGS2phylo/metadata_input_Campylobacter.txt --with-fasta --output-prefix metadata_out/metadata_required_Campylobacter

python3 bin/WGS2Phylo.py --output_dir data/example_data/WGS2phylo/results_Campylobacter/ --organism campylobacter --supplemental-file data/example_data/WGS2phylo/metadata_input_Campylobacter.txt --with-fasta --extra-fields --output-prefix metadata_out/metadata_expanded_Campylobacter

```

### RSV

```
python3 bin/WGS2Phylo.py --output_dir data/example_data/WGS2phylo/results_rsv/ --organism rsv --supplemental-file data/example_data/WGS2phylo/metadata_input_rsv.txt --with-fasta --output-prefix metadata_out/metadata_required_rsv

python3 bin/WGS2Phylo.py --output_dir data/example_data/WGS2phylo/results_rsv/ --organism rsv --supplemental-file data/example_data/WGS2phylo/metadata_input_rsv.txt --with-fasta --extra-fields --output-prefix metadata_out/metadata_expanded_rsv
```

## Tests
Go to tests/WGS2Phylo and execute:
```
pytest test_WGS2Phylo.py --data-dir ../../data/example_data/WGS2phylo/unit_tests/ -v
``` 

-----------

# MST Tree

The bacterial pipeline produces an interactive **Minimum Spanning Tree (MST)** in HTML format for all samples listed in the metadata file.
The MST is constructed based on allelic differences observed between profiles, adjusted for missing loci (following the pHierCC methodology).

Alongside the HTML plot the pipeline writes:
- `*_MST.tsv` – the MST edge list (`source`, `target`, `distance` in allelic differences),
- `*_MST.nwk` – a sample-level Newick representation of the same MST, added as a third tree tab (`cgMLST MST`) in the Microreact project, next to the phylogenetic tree and the time tree.

All trees share one Microreact pane and are switched with the tabs above it. Viral projects, which have no MST, show the same pane with two tabs.

The MST itself is calculated between unique cgMLST sequence types (ST). Because Microreact links tree tips to metadata rows by sample identifier, every ST becomes an internal node in the Newick file and its samples are attached as zero-length terminal branches. Branch lengths between ST nodes are the allelic distances.

The visualization root is the deterministic weighted graph centre: the ST that minimizes the largest allelic distance to any other ST. **This is a visualization root only and does not represent inferred ancestry.**

Samples whose cgMLST profile cannot be resolved against the profile database are reported in the process log and omitted from the Newick file; their metadata rows are still available in the Microreact map and table.

As a result, the bacterial pipeline shell wrapper needs a `--db` directory providing information about alleles identified at each locus for a given Sequence Type (ST) in the cgMLST schema. It defaults to the shared `/mnt/unity_nfs/external_databases`.

## External Database Structure

The external database directory must contain cgMLST profile definitions for supported species.  
Each species directory includes:
- a main `profiles.list` file (publicly released cgMLST profiles)
- a `local/profiles_local.list` file (internal or temporary profiles, optional)

### Expected layout

External datases strucutre is predefined and described in details in our [Sequnecing pipline](https://github.com/mkadlof/pzh_pipeline_viral), however the only requiered
files are one with profiles information for supported species. Follwoing structure of the `--db` directory (by default `/mnt/unity_nfs/external_databases`) must be respected

```

/mnt/unity_nfs/external_databases/
├── cgmlst/
│   ├── Salmonella/
│   │   ├── profiles.list
│   │   └── local/
│   │       └── profiles_local.list
│   │
│   ├── Escherichia/
│   │   ├── profiles.list
│   │   └── local/
│   │       └── profiles_local.list
│   │
│   └── Campylobacter/
│       └── jejuni/
│           ├── profiles.list
│           └── local/
│               └── profiles_local.list

```

`profiles_local.list` includes additional profiles that might not be present in "main" file (e.g. temporal identifier, private STs etc.)

## Execution

Plotting of the MST is integrated*into the bacterial pipeline.No additional steps are required.  
The resulting HTML file will be saved in the pipeline’s output directory.

## Tests

Run from the repository root:
```bash
pytest tests/MST_bacteria -v
```

Dependencies: numpy, pandas, scipy, networkx, Biopython, Plotly, Click, and pytest

