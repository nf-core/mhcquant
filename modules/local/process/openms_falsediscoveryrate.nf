// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process OPENMS_FALSEDISCOVERYRATE  {

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(id_file_idx)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_idx_fdr.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    when:
        !params.skip_quantification

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