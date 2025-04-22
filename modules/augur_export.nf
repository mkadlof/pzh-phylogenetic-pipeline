process augur_export {
    input:
    path tree
    path metadata
    path node_data

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