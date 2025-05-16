process treetime {
    input:
    tuple val(segmentId), path(alignment), path(tree)
    path metadata

    output:
    tuple val(segmentId), path("timetree.nwk"), path("*.node_data.json")

    script:
    """
    augur refine \
        --alignment ${alignment} \
        --tree ${tree} \
        --metadata ${metadata} \
        --output-tree timetree.nwk \
        --keep-polytomies \
        --branch-length-inference joint \
        --keep-root \
        --timetree \
        --verbosity 6
    """
}