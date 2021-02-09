// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process RUN_PERCOLATOR {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(id_file_psm)
        val fdr_level

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_all_ids_merged_psm_perc.idXML")

    if (params.klammer && params.description_correct_features == 0) {
        log.warn('Klammer was specified, but description of correct features was still 0. Please provide a description of correct features greater than 0.')
        log.warn('Klammer has been turned off!')
    }

    script:
        if (params.description_correct_features > 0 && params.klammer){
            """
                OMP_NUM_THREADS=${task.cpus} \\
                PercolatorAdapter -in ${id_file_psm} \\
                    -out ${Sample}_all_ids_merged_psm_perc.idXML \\
                    -seed 4711 \\
                    -trainFDR 0.05 \\
                    -testFDR 0.05 \\
                    -enzyme no_enzyme \\
                    $fdr_level \\
                    -subset-max-train ${params.subset_max_train} \\
                    -doc ${params.description_correct_features} \\
                    -klammer
            """
        } else {
            """
                OMP_NUM_THREADS=${task.cpus} \\
                PercolatorAdapter -in ${id_file_psm} \\
                    -out ${Sample}_all_ids_merged_psm_perc.idXML \\
                    -seed 4711 \\
                    -trainFDR 0.05 \\
                    -testFDR 0.05 \\
                    -threads ${task.cpus} \\
                    -enzyme no_enzyme \\
                    $fdr_level \\
                    -subset-max-train ${params.subset_max_train} \\
                    -doc ${params.description_correct_features} \\
            """
        }
}