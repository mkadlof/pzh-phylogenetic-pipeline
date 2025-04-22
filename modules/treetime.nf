process treetime {
    input:
    path alignment
    path metadata
    path tree

    output:
    path "timetree.nwk", emit: timetree
    path "*.node_data.json", emit: node_data

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