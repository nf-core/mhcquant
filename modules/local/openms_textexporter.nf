// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_TEXTEXPORTER {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0"
    } else {
        container "quay.io/biocontainers/openms:2.6.0--h4afb90d_0"
    }

    input:
        tuple val(meta), path(consensus_resolved)

    output:
        tuple val(meta), path("*.tsv"), emit: tsv
        path "versions.yml"           , emit: versions

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${meta.id}_${options.suffix}" : "${meta.id}"

        """
        TextExporter -in $consensus_resolved \\
            -out ${prefix}.tsv \\
            -threads $task.cpus \\
            -id:add_hit_metavalues 0 \\
            -id:add_metavalues 0 \\
            -id:peptides_only

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
