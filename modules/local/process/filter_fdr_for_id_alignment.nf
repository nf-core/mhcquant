// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process FILTER_FDR_FOR_ID_ALIGNMENT {

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(id_file_idx_fdr) 

    output:
        tuple val(id), val("$Sample"), val(Condition), file("${id}_-_${Sample}_${Condition}_idx_fdr_filtered.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    when:
        !params.skip_quantification

    script:
    """
        IDFilter -in ${id_file_idx_fdr} \\
            -out ${id}_-_${Sample}_${Condition}_idx_fdr_filtered.idXML \\
            -threads ${task.cpus} \\
            -score:pep ${params.fdr_threshold} \\
            -precursor:length '${params.peptide_min_length}:${params.peptide_max_length}' \\
            -remove_decoys \\
            -delete_unreferenced_peptide_hits

        FileInfo --help &> openms.version.txt
    """
}