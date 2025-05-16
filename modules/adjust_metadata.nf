process adjust_metadata {
    input:
    path metadata

    output:
    path "influenza_metadata_adjusted.tsv", emit: adjusted_metadata

    script:
    """
    adjust_metadata.py ${metadata}
    """
}