process augur_filter_sequences {
    input:
    path fasta
    path index_csv
    path metadata
    path exclude

    output:
    path "valid_strains.fasta", emit: out

    script:
    """
    augur filter \
        --sequences ${fasta} \
        --sequence-index ${index_csv} \
        --metadata ${metadata} \
        --exclude ${exclude} \
        --output valid_strains.fasta
    """
}