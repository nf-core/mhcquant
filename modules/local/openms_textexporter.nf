// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process OPENMS_TEXTEXPORTER {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(Sample), path(consensus_resolved) 

    output:
        tuple val(Sample), path("${Sample}.csv"), emit: csv   
        path  "*.version.txt", emit: version 

    script:
        def software = getSoftwareName(task.process)

        """
            TextExporter -in ${consensus_resolved} \\
                -out ${Sample}.csv \\
                -threads ${task.cpus} \\
                -id:add_hit_metavalues 0 \\
                -id:add_metavalues 0 \\
                -id:peptides_only
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}