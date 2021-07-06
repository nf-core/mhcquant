// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process OPENMS_FALSEDISCOVERYRATE  {

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), path(id_file_idx)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), path("${Sample}_${Condition}_${id}_idx_fdr.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)

        """
            FalseDiscoveryRate -in ${id_file_idx} \\
                -protein 'false' \\
                -out ${Sample}_${Condition}_${id}_idx_fdr.idXML \\
                -threads ${task.cpus}
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}