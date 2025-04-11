process augur_index_sequences {
    input:
    path fasta

    output:
    path "sequence-index.csv", emit: out

    script:
    """
    augur index --sequences ${fasta} --output sequence-index.csv
    """
}