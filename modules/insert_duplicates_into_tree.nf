process insert_duplicates_into_tree {
    input:
    path tree
    path ids

    output:
    path "consensus_tree.nwk"

    script:
    """
    insert_missing_duplicated_sequences_into_tree.py \
        --tree ${tree} \
        --ids ${ids} > consensus_tree.nwk
    """
}