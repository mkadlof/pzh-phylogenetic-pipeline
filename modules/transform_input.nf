process transform_input {
    input:
    path input_dir

    script:
    """
    # Read segments in the first sample directory
    segments=(\$(ls -1 "\$(ls -d ${input_dir}/*/ | head -n1)" ))

    samples=(\$(ls -d ${input_dir}/*/))

    for segment in \${segments[@]}; do
        # Create empty fasta for segment
        touch \${segment}

        # Concatenate all samples for this segment
        for sample in \${samples[@]}; do
            cat "\${sample}\${segment}" >> \${segment}
        done
    done
    """
}