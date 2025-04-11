// Variables for processes
def hostname = "hostname".execute().text.trim() // We need that to overwrite a default "container" options from config, used by the alphafold
ExecutionDir = new File('.').absolutePath

// ALL parameters are setup usomg bash wrapper except enterobase_api_token that MUST be part of nextflow config
// Comments were preserved in the  nf file for a local executor
params.input_dir = ""
params.input_type = ""
params.main_image = "" 
params.results_dir = ""
params.prokka_image = ""
params.threads = 1
params.metadata = ""


// User must use our config that has two profiles slurm and local, nextflow must be initialized with one of them

if ( !workflow.profile || ( workflow.profile != "slurm" && workflow.profile != "local") ) {
   println("Nextflow run must be executed with -profile option. The specified profile must be either \"local\" or \"slurm\".")
   System.exit(1)
}

// QC params
params.threshold_Ns = 100
params.threshold_ambiguities = 100


// Processes 


process run_prokka {
  // publishDir "${params.results_dir}/${x}/", mode: 'copy', pattern: "${x}_prokka*"
  container  = params.prokka_image
  tag "Predicting genes for sample $x"
  cpus { params.threads > 25 ? 25 : params.threads }
  memory "10 GB"
  time "20m"
  input:
  tuple val(x), path(fasta)
  output:
  path("${x}_prokka.gff")
  script:
  """
  prokka --metagenome --cpus ${task.cpus} --outdir prokka_out --prefix prokka_out --compliant --kingdom Bacteria $fasta
  echo -e "{\\"status\\": \\"tak\\", \
            \\"prokka_gff\\": \\"${params.results_dir}/${x}/${x}_prokka.gff\\", \
            \\"prokka_ffn\\": \\"${params.results_dir}/${x}/${x}_prokka.ffn\\"}" >> prokka.json
  # Following files are usefull for phylogenetic analyis
  mv prokka_out/prokka_out.gff ${x}_prokka.gff
  mv prokka_out/prokka_out.ffn ${x}_prokka.ffn

  """
}

process run_roary {
  container  = params.main_image
  tag "Predicting pangenome with roary"
  cpus { params.threads > 25 ? 25 : params.threads }
  memory "30 GB"
  time "1h"
  input:
  path(gff)
  output:
  tuple path("core_genes_alignment.fasta"), path("core_genes_alignment.embl")
  script:
  """
  # -f to nazwa katalogu z outputem
  # -e create a multiFASTA alignment
  # -n fast core gene alignment with MAFFT, use with -e
  # -v verbose
  # -p liczba rdzeni
  # -i minimum percentage identity for blastp [95]
  # -cd FLOAT percentage of isolates a gene must be in to be core [99]
  roary -p ${task.cpus} -i 95 -cd 95  -f ./roary_output -e -n *.gff
  cp roary_output/core_gene_alignment.aln core_genes_alignment.fasta
  cp roary_output/core_alignment_header.embl core_genes_alignment.embl
  """
}

process augur_index_sequences {
    container  = params.main_image
    tag "Indexing sequences with augur"
    cpus 1
    memory "30 GB"
    time "1h"
    
    input:
    tuple path(fasta), path(embl)

    output:
    tuple path(fasta), path(embl), path("index.csv")

    script:
    """
    augur index --sequences ${fasta} --output index.csv
    """
}

process augur_filter_sequences {
    container  = params.main_image
    tag "Filtering out sequences with augur"
    cpus 1
    memory "30 GB"
    time "1h"

    input:
    tuple path(fasta), path(embl), path(index)
    path(metadata)
    output:
    tuple path("valid_sequences.fasta"), path(embl), path(metadata)

    script:
    """
    # For NOW we are liberal when it comes to sequences quality
    # the script only checks columns 5 and 6 in $index i.e. Ns and abigous
    python /opt/docker/custom_scripts/filter_low_quality_sequences.py --output_dir . \
                                                                      --threshold_Ns ${params.threshold_Ns} \
                                                                      --threshold_ambiguities ${params.threshold_ambiguities} \
                                                                      $index
    
    # For now we use augur filter to preprare fasta file without invalid_strains.txt prepared with filter_low_quality_sequences script 
    # Other usefull options --min-length --max-length  --group-by which we do not use for now
    augur filter \
        --sequences ${fasta} \
        --sequence-index ${index} \
        --metadata ${metadata} \
        --exclude invalid_strains.txt \
        --output-sequences valid_sequences.fasta

    """
}

process save_input_to_log {
  tag "Dummy process"
  cpus 1
  memory "1 GB"
  time "1m"
  input:
  path(x)
  output:
  stdout
  script:
  """
  echo ${x} >> log
  """
}

// MAIN WORKFLOW //

workflow {
Channel
    .fromPath("${params.metadata}")
    .set {metadata_channel}

// Prepare gff input

if (params.input_type == 'fasta') {
    Channel
        .fromPath("${params.input_dir}/*")
        .map { file -> tuple(file.getName().split("\\.")[0], file) }
        .set { initial_fasta }

    gff_input = run_prokka(initial_fasta).collect()

} else if (params.input_type == 'gff') {
    Channel
        .fromPath("${params.input_dir}/*")
        .collect()
        .set { gff_input }

} else {
    println("--input_type must be either fasta or gff")
    System.exit(1)
}

roary_out = run_roary(gff_input)

augur_index_sequences_out = augur_index_sequences(roary_out)

augur_filter_sequences_out = augur_filter_sequences(augur_index_sequences_out, metadata_channel)

// save_input_to_log(gff_input)

}
