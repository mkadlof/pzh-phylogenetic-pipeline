process insert_duplicates_into_alignment {
    input:
    path alignment
    path ids

    output:
    path "alignment_with_duplicates.fasta", emit: out

    script:
    """
    add_or_remove_duplicates_from_alignment.py \
        --alignment ${alignment} \
        --duplicated_ids ${ids} \
        --action insert \
        --output ./alignment_with_duplicates.fasta
    """
}