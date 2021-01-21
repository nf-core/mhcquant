// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process FILTER_FDR_FOR_ID_ALIGNMENT {

    input:
    tuple val(id), val(Sample), val(Condition), file(id_file_idx_fdr) 

    output:
    tuple val(id), val("$Sample"), val(Condition), file("${id}_-_${Sample}_${Condition}_idx_fdr_filtered.idXML")

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
     """
// 
}