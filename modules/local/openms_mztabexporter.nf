// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process OPENMS_MZTABEXPORTER {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }
    
    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), path(mztab)

    output:
        tuple val("$id"), val("$Sample"), path("*.mzTab"), emit: mztab   
        path  "*.version.txt", emit: version

    script:
        def prefix = options.suffix ? "${Sample}_${options.suffix}" : "${Sample}"
        def software = getSoftwareName(task.process)
        
        """
            MzTabExporter -in ${mztab} \\
                -out ${prefix}.mzTab \\
                -threads ${task.cpus}
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}