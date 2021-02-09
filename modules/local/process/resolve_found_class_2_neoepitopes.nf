// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process RESOLVE_FOUND_CLASS_2_NEOEPITOPES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=2.7.15" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:2.7.15"
    } else {
        container "quay.io/biocontainers/python:2.7.15"
    }

    echo true

    input:
        tuple val(id), val(Sample), file(mztab), file(neoepitopes)

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_found_neoepitopes_class_2.csv")

    when:
        params.include_proteins_from_vcf & params.predict_class_2

    script:
        """
        resolve_neoepitopes.py -n ${neoepitopes} -m ${mztab} -f csv -o ${Sample}_found_neoepitopes_class_2
        """
}