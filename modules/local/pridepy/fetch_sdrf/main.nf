process PRIDEPY_FETCH_SDRF {
    label 'process_single'
    tag "$pride_id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pridepy:0.0.12--pyhdfd78af_0' :
        'quay.io/biocontainers/pridepy:0.0.12--pyhdfd78af_0' }"

    input:
    val pride_id

    output:
    path "*.sdrf.tsv" , emit: sdrf
    tuple val("${task.process}"), val('pridepy'), eval("pip show pridepy 2>/dev/null | grep Version | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pridepy stream-files-metadata -a "${pride_id}" -o files_metadata.json
    sdrf_name=\$(python3 -c "import json,sys; sdrfs=[f['fileName'] for f in json.load(open('files_metadata.json')) if f['fileName'].endswith('.sdrf.tsv')]; print(sdrfs[0]) if sdrfs else sys.exit('ERROR: No SDRF file found for ${pride_id}')")
    pridepy download-file-by-name -a "${pride_id}" -f "\$sdrf_name" -o . -p ftp
    """

    stub:
    """
    touch ${pride_id}.sdrf.tsv
    """
}
