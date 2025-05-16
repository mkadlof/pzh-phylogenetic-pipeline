process augur_align {
    tag "${segmentId}"
    input:
    tuple val(segmentId), path(fasta)

    output:
    tuple val(segmentId), path("aligned.fasta"), emit: out

    script:
    """
    augur align --sequences ${fasta} --output aligned.fasta
    """
}