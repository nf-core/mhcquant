// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_IDFILTER {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0"
    } else {
        container "quay.io/biocontainers/openms:2.6.0--h4afb90d_0"
    }

    input:
        tuple val(meta), path(idxml), file(peptide_filter)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml", emit: versions

    script:
        def software = getSoftwareName(task.process)
        def whitelist = "${peptide_filter}"
        def prefix = options.suffix ? "${idxml.baseName}_${options.suffix}" : "${meta.id}_-_idx_fdr_filtered"

        if (whitelist == "input.2") {
            whitelist = " "
        }

        """
        IDFilter -in $idxml \\
            -out ${prefix}.idXML \\
            -threads $task.cpus \\
            $options.args \\
            $whitelist
        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
