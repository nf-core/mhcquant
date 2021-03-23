// Import generic module functions
include { initOptions } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process ALIGN_IDS {
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(id_names)

    output:
        tuple val("$Sample"), file("*.trafoXML"), emit: trafoxml   
        path  "*.version.txt", emit: version

    when:
        !params.skip_quantification

    script:
        def out_names = id_names.collect { it.baseName+'.trafoXML' }.join(' ')

        """
            MapAlignerIdentification 
            -in $id_names \\
            -trafo_out $out_names \\
            -model:type linear \\
            -algorithm:max_rt_shift ${params.max_rt_alignment_shift}

            FileInfo --help &> openms.version.txt
        """
}