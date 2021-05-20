// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

//TODO: combine in a subflow --> when needs to be removed
process OPENMS_MAPRTTRANSFORMER {

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(file_align), file(file_trafo)

    output:
            
        tuple val("$id"), val("$Sample"), val("$Condition"), file("*_aligned.*"), emit: aligned   
        path  "*.version.txt", emit: version
            
    when:
        !params.skip_quantification

    script:
        def fileExt = file_align.collect { it.name.tokenize("\\.")[1] }.join(' ')
        def software = getSoftwareName(task.process)

        """
            MapRTTransformer -in ${file_align} \\
                -trafo_in ${file_trafo} \\
                -out ${Sample}_${Condition}_${id}_aligned.${fileExt} \\
                -threads $task.cpus
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}