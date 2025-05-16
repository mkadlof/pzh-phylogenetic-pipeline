process augur_filter_sequences {
    input:
    tuple val(segmentId), path(fasta), path(index_csv), path(exclude)
    path metadata

    output:
    tuple val(segmentId), path("valid_strains.fasta")

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