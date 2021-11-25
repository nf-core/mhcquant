// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_PEPTIDEINDEXER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0"
    } else {
        container "quay.io/biocontainers/openms:2.6.0--h4afb90d_0"
    }

    input:
        tuple val(meta), path(idxml), path(fasta)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml", emit: versions

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${idxml.baseName}_${options.suffix}" : "${idxml.baseName}_idx"

        """
            PeptideIndexer -in ${idxml} \\
                    -out ${prefix}.idXML \\
                    -threads ${task.cpus} \\
                    -fasta ${fasta} \\
                    -decoy_string DECOY \\
                    -enzyme:specificity none

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
            END_VERSIONS
        """
}
