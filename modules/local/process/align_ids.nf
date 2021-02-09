// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process ALIGN_IDS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(id_names)

    output:
        tuple val("$Sample"), file("*.trafoXML")

    when:
        !params.skip_quantification

    script:
        def out_names = id_names.collect { it.baseName+'.trafoXML' }.join(' ')

        """
            MapAlignerIdentification -in $id_names \\
            -trafo_out $out_names \\
            -model:type linear \\
            -algorithm:max_rt_shift ${params.max_rt_alignment_shift}
        """
}