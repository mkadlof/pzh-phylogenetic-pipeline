process remove_duplicates_from_alignment {
    input:
    path alignment
    path duplicated_ids

    output:
    path "deduplicated_alignment.fasta", emit: deduplicated_alignment

    script:
    """
    remove_duplicates_from_alignment.py -a ${alignment} -i ${duplicated_ids} -o deduplicated_alignment.fasta
    """
}