process transform_input {
    cpus 1
    memory "30 GB"
    time "1h"
    input:
    path input_dir

    output:
    path "*.fasta", emit: fastas

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