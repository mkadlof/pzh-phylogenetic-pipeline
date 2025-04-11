process find_identical_sequences {
    input:
    path fasta

    output:
    path "valid_strains_unique.fasta", emit: uniq_fasta
    path "valid_strains_ident_seq.csv", emit: duplicated_ids

    script:
    """
    find_identical_sequences.py --input ${fasta} --output_dir .
    """
}