process augur_index_sequences {
    container = params.main_image
    tag "Concatanating data"
    cpus 1
    memory "30 GB"
    time "1h"
    input:
    input:
    path fasta
    
    output:
    path "sequence-index.csv", emit: sequence_index

    script:
    """
    augur index --sequences ${fasta} --output sequence-index.csv
    """
}
