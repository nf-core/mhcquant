// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process CALCULATE_FDR_FOR_ID_ALIGNMENT  {
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
        tuple val(id), val(Sample), val(Condition), file(id_file_idx)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_idx_fdr.idXML") 

    when:
        !params.skip_quantification

    script:
        """
            FalseDiscoveryRate -in ${id_file_idx} \\
            -protein 'false' \\
            -out ${Sample}_${Condition}_${id}_idx_fdr.idXML \\
            -threads ${task.cpus}
        """
}