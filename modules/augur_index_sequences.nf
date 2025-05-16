process augur_index_sequences {
    tag "${segmentId}"
    input:
    tuple val(segmentId), path(fasta)

    output:
    tuple val(segmentId), path("sequence-index.csv"), emit: sequence_index

    script:
    """
    augur index --sequences ${fasta} --output sequence-index.csv
    """
}