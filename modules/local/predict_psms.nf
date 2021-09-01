// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREDICT_PSMS {
    tag "$meta"
    label 'process_medium'

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
        tuple val(meta), path(perc_mztab), path(psm_mztab), val(allotypes)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path  "*.version.txt", emit: version

    script:
        def prefix = options.suffix ? "${meta.id}_${options.suffix}" : "${meta.id}_peptide_filter"

    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab_for_filtering.py ${params.subset_affinity_threshold} '${allotypes}' ${perc_mztab} ${psm_mztab} ${prefix}.idXML
        mhcflurry-predict --version &> mhcflurry.version.txt
    """

}
