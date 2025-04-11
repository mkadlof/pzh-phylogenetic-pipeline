// Variables for processes
def hostname = "hostname".execute().text.trim() // We need that to overwrite a default "container" options from config, used by the alphafold
ExecutionDir = new File('.').absolutePath

// ALL parameters are setup usomg bash wrapper except enterobase_api_token that MUST be part of nextflow config
// Comments were preserved in the  nf file for a local executor
params.input_dir = ""
params.input_type = ""
params.some_parameter = ""
params.main_image = "" 
params.results_dir = ""
params.prokka_image = ""

// User must use our config that has two profiles slurm and local, nextflow must be initialized with one of them

if ( !workflow.profile || ( workflow.profile != "slurm" && workflow.profile != "local") ) {
   println("Nextflow run must be executed with -profile option. The specified profile must be either \"local\" or \"slurm\".")
   System.exit(1)
}

// Processes 


process run_prokka {
  // publishDir "${params.results_dir}/${x}/", mode: 'copy', pattern: "${x}_prokka*"
  container  = params.prokka_image
  tag "Predicting genes for sample $x"
  cpus { params.threads > 25 ? 25 : params.threads }
  memory "10 GB"
  time "20m"
  input:
  tuple val(x), path(fasta), val(QC_status)
  output:
  path("${x}_prokka.gff")
  script:
  """
  if [[ ${QC_status} == "nie"  ]]; then
    mkdir prokka_out; touch prokka_out/prokka_out_dummy.gff; touch prokka_out/prokka_out_dummy.faa; touch prokka_out/prokka_out_dummy.ffn; touch prokka_out/prokka_out_dummy.tsv
    ERROR_MSG="Initial QC received by this module was nie"
    echo -e "{\\"status\\": \\"nie\\", \
              \\"error_message\\": \\"\${ERROR_MSG}\\"}"  >> prokka.json
    # json z informacja o bledzie jakosci
  else
      prokka --metagenome --cpus ${task.cpus} --outdir prokka_out --prefix prokka_out --compliant --kingdom Bacteria $fasta
      echo -e "{\\"status\\": \\"tak\\", \
            \\"prokka_gff\\": \\"${params.results_dir}/${x}/${x}_prokka.gff\\", \
            \\"prokka_ffn\\": \\"${params.results_dir}/${x}/${x}_prokka.ffn\\"}" >> prokka.json
  fi
  # Following files are usefull for phylogenetic analyis
  mv prokka_out/prokka_out.gff ${x}_prokka.gff
  mv prokka_out/prokka_out.ffn ${x}_prokka.ffn

  """
}

process save_input_to_log {
  tag "Dummy process"
  cpus 1
  memory "1 GB"
  time "1m"
  input:
  tuple val(x)
  output:
  stdout
  script:
  """
  echo ${x} >> log
  """
}

// MAIN WORKFLOW //

workflow {

// Prepare gff input
if(params.input_type == 'fasta') {
Channel
  .fromPath(params.input_dir/*)
  .map {it -> tuple(it.getName().split("\\.")[0], it)}
  .set {initial_fasta}

gff_input =  run_prokka(initial_fasta).collect()

} elif(params.input_type == 'gff') {
Channel
  .fromPath(params.input_dir/*)
  .collect()
  .set {gff_input}


} else {
  println("--input_type must be either fasta or gff")
  System.exit(1)
}
   

save_input_to_log(gff_input)
// run roary to get common genes 
// run_roary_out = run_roary(gff_input)

}
