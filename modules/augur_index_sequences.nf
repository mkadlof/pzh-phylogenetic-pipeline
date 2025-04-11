process augur_index_sequences {
    input:
    path fasta

    output:
    path "sequence-index.csv", emit: sequence_index

    script:
    """
    augur index --sequences ${fasta} --output sequence-index.csv
    """
}