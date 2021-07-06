// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process OPENMS_IDMERGER {
    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), path(ids_aligned)

    output:
        tuple val("$id"), val("$Sample"), val(Condition), path("*.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)

        """
            IDMerger -in $ids_aligned \\
                -out ${Sample}_all_ids_merged.idXML \\
                -threads ${task.cpus} \\
                -annotate_file_origin \\
                -merge_proteins_add_PSMs
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}