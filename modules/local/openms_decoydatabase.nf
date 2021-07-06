// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process OPENMS_DECOYDATABASE {
    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), path(fastafile)

    output:
        tuple val("$id"), val("$Sample"), path("${fastafile.baseName}_decoy.fasta"), emit: decoy
        path "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)

        """
            DecoyDatabase -in ${fastafile} \\
                -out ${fastafile.baseName}_decoy.fasta \\
                -decoy_string DECOY_ \\
                -decoy_string_position prefix
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}