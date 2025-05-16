process insert_duplicates_into_alignment {
    tag "${segmentId}"
    input:
    tuple val(segmentId), path(alignment), path(ids)

    output:
    tuple val(segmentId), path("alignment_with_duplicates.fasta")

    script:
    """
    add_or_remove_duplicates_from_alignment.py \
        --alignment ${alignment} \
        --duplicated_ids ${ids} \
        --action insert \
        --output ./alignment_with_duplicates.fasta
    """
}