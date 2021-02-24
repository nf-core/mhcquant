// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process ALIGN_MZML_FILES {
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(mzml_file_align), file(id_file_trafo)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_aligned.mzML"), emit: mzml   
        path  "*.version.txt", emit: version

    when:
        !params.skip_quantification

    script:
        """
            MapRTTransformer -in ${mzml_file_align} \\
            -trafo_in ${id_file_trafo} \\
            -out ${Sample}_${Condition}_${id}_aligned.mzML \\
            -threads ${task.cpus}
            
            FileInfo --help &> openms.version.txt
        """
}