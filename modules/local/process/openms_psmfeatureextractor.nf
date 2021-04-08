// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process OPENMS_PSMFEATUREEXTRACTOR {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(id_file_merged)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_all_ids_merged_psm.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        """
            PSMFeatureExtractor -in ${id_file_merged} \\
            -out ${Sample}_all_ids_merged_psm.idXML \\
            -threads ${task.cpus} 

            FileInfo --help &> openms.version.txt
        """
}