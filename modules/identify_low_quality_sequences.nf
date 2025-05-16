process identify_low_quality_sequences {
    tag "${segmentId}"
    cpus 1
    memory "30 GB"
    time "1h"
    input:
    tuple val(segmentId), path(index_csv)

    output:
    tuple val(segmentId), path("invalid_strains.txt"), emit: out

    script:
    """
    identify_low_quality_sequences.py \
        --output_dir . \
        --threshold_Ns ${params.threshold_Ns} \
        --threshold_ambiguities ${params.threshold_ambiguities} \
        ${index_csv}
    """
}