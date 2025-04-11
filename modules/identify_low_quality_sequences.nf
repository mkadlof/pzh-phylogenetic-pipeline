process identify_low_quality_sequences {
    input:
    path index_csv

    output:
    path "invalid_strains.txt", emit: out

    script:
    """
    pwd
    echo $PATH
    identify_low_quality_sequences.py --output_dir . ${index_csv}
    """
}