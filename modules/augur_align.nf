process augur_align {
    input:
    path fasta

    output:
    path "aligned.fasta", emit: out

    script:
    """
    augur align --sequences ${fasta} --output aligned.fasta
    """
}