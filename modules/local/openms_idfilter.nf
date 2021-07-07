// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

options = initOptions(params.options)

process OPENMS_IDFILTER {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), path(id_file), file(peptide_filter)

    output:
        tuple val(id), val(Sample), val(Condition), path("*.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${Sample}_${options.suffix}" : "${id}_-_${Sample}_${Condition}_idx_fdr_filtered"
        def whitelist = "${peptide_filter}"

        if (whitelist == "input.2") {
            whitelist = " "
        } 

        """
            IDFilter -in ${id_file} \\
                -out ${prefix}.idXML \\
                -threads ${task.cpus} \\
                $options.args ${whitelist}
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}