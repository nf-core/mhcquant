process PRIDE_DOWNLOAD_FILE {
    label 'process_single'
    tag "${file_name}"

    conda "bioconda::pridepy=0.0.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pridepy:0.0.12--pyhdfd78af_0' :
        'quay.io/biocontainers/pridepy:0.0.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(file_name), val(pride_accession)

    output:
    tuple val(meta), path("${file_name}"), emit: downloaded_file
    tuple val("${task.process}"), val('pridepy'), eval("pip show pridepy 2>/dev/null | grep Version | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pridepy download-file-by-name \\
        -a "${pride_accession}" \\
        -f "${file_name}" \\
        -o . \\
        -p ftp

    # pridepy exits 0 even on download failure — validate file is non-empty
    [ -s "${file_name}" ] || { echo "ERROR: Downloaded file ${file_name} is empty"; exit 1; }
    """

    stub:
    """
    touch "${file_name}"
    """
}
