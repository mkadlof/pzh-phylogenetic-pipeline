process adjust_metadata {
    input:
    path metadata

    output:
    path "adjusted_metadata.tsv", emit: adjusted_metadata

    script:
    """
    # Module to be implementd.
    cat ${metadata} > adjusted_metadata.tsv
    """
}