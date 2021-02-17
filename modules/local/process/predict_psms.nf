// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process PREDICT_PSMS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::mhcflurry=1.4.3--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcflurry:1.4.3--py_0"
    } else {
        container "quay.io/biocontainers/mhcflurry:1.4.3--py_0"
    }

    input:
        tuple val(Sample), val(id), file(perc_mztab_file), file(psm_mztab_file), val(d), val(allotypes_refine)

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_peptide_filter.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    when:
        params.refine_fdr_on_predicted_subset

    script:
    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab_for_filtering.py ${params.subset_affinity_threshold} '${allotypes_refine}' ${perc_mztab_file} ${psm_mztab_file} ${Sample}_peptide_filter.idXML

        mhcflurry-predict --version &> mhcflurry.version.txt
    """
}
