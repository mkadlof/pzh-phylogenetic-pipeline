process augur_export {
    tag "${segmentId}"
    input:
    tuple val(segmentId), path(tree), path(node_data)
    path metadata

    output:
    path "auspice.json", emit: out

    script:
    """
    augur export v2 \
        --tree ${tree} \
        --metadata ${metadata} \
        --auspice-config /etc/auspice/auspice_config.json \
        --title "PZH viral phylogenetic pipeline" \
        --node-data ${node_data} \
        --output auspice.json
    """
}