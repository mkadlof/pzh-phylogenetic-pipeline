process find_identical_sequences {
    input:
    path fasta

    output:
    path "valid_strains_unique.fasta", emit: out
    path "valid_strains_ident_seq.csv", emit: ids

    script:
    """
    find_identical_sequences.py --input ${fasta} --output_dir .
    """
}