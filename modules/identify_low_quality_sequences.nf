process identify_low_quality_sequences {
    input:
    tuple val(segmentId), path(index_csv)

    output:
    tuple val(segmentId), path("invalid_strains.txt"), emit: out

    script:
    """
    identify_low_quality_sequences.py --output_dir . ${index_csv}
    """
}