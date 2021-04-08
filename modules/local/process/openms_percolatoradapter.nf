// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.option = [:]
options = initOptions(params.options)

//TODO: combine in a subflow --> when needs to be removed
process RUN_PERCOLATOR_ON_PREDICTED_SUBSET {
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
        tuple val(id), val(Sample), file(id_file_psm_subset)
        val fdr_level

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_perc_subset.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    when:
        params.refine_fdr_on_predicted_subset

    script:
        """
        OMP_NUM_THREADS=${task.cpus} \\
        PercolatorAdapter -in ${id_file_psm_subset} \\
            -out ${Sample}_perc_subset.idXML \\
            -seed 4711 \\
            -trainFDR 0.05 \\
            -testFDR 0.05 \\
            -enzyme no_enzyme \\
            -subset-max-train ${params.subset_max_train} \\
            -doc ${params.description_correct_features} \\
            $fdr_level

        FileInfo --help &> openms.version.txt
        """
}