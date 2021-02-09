// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process MERGE_ALIGNED_IDMXL_FILES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(ids_aligned)

    output:
        tuple val("$id"), val("$Sample"), val(Condition), file("${Sample}_all_ids_merged.idXML")

    script:
    """
        IDMerger -in $ids_aligned \\
        -out ${Sample}_all_ids_merged.idXML \\
        -threads ${task.cpus}  \\
        -annotate_file_origin  \\
        -merge_proteins_add_PSMs
    """
}