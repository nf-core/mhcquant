// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process PREDICT_PSMS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(Sample), val(id), path(perc_mztab_file), path(psm_mztab_file), val(d), val(allotypes_refine)

    output:
        tuple val("$id"), val("$Sample"), path("${Sample}_peptide_filter.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab_for_filtering.py ${params.subset_affinity_threshold} '${allotypes_refine}' ${perc_mztab_file} ${psm_mztab_file} ${Sample}_peptide_filter.idXML
        mhcflurry-predict --version &> mhcflurry.version.txt
    """

}
