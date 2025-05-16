process augur_align {
    tag "${segmentId}"
    cpus 1
    memory "30 GB"
    time "1h"
    input:
    tuple val(segmentId), path(fasta)

    output:
    tuple val(segmentId), path("aligned.fasta"), emit: out

    script:
    """
    augur align --sequences ${fasta} --output aligned.fasta
    """
}