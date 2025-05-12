process transform_input {
    container = params.main_image
    tag "Concatanating data"
    cpus 1
    memory "30 GB"
    time "1h"
    input:
    path(x) 
    output:
    tuple path("HA.fasta"), path("NA.fasta"), path("PB1.fasta")
    script:
    """
    # Extract fasta file from input fasta 
    python3 /opt/docker/extract_segment.py --input_dir . --segment_name HA
    python3 /opt/docker/extract_segment.py --input_dir . --segment_name NA
    python3 /opt/docker/extract_segment.py --input_dir . --segment_name PB1
    # and so on ...
    """
}


