process insert_duplicates_into_tree {
    tag "${segmentId}"
    input:
    tuple val(segmentId), path(tree), path(ids)

    output:
    tuple val(segmentId), path("consensus_tree.nwk")

    script:
    """
    insert_missing_duplicated_sequences_into_tree.py \
        --tree ${tree} \
        --ids ${ids} > consensus_tree.nwk
    """
}