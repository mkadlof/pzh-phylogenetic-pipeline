process transform_input {
    input:
    path input_dir

    output:
    path "*.fasta", emit: fastas
//    tuple val("PB2"), path("PB2.fasta")
//    tuple val("PB1"), path("PB1.fasta")
//    tuple val("PA"),  path("PA.fasta")
//    tuple val("HA"),  path("HA.fasta")
//    tuple val("NP"),  path("NP.fasta")
//    tuple val("NA"),  path("NA.fasta")
//    tuple val("MP"),  path("MP.fasta")
//    tuple val("NS"),  path("NS.fasta")

    script:
    """
    samples=(\$(ls -d ${input_dir}/*/))
    segments=(PB2 PB1 PA HA NP NA MP NS)

    for idx in \$(seq 1 8); do
        segment=\${segments[\$((idx-1))]}
        infile="output_chr\${idx}_\${segment}.fasta"
        outname="\${segment}.fasta"

        for sample in "\${samples[@]}"; do
            cat "\${sample}/\${infile}" >> \$outname
        done
    done
    """
}