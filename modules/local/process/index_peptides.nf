// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process INDEX_PEPTIDES {

    input:
    tuple val(Sample), val(id), val(Condition), file(id_file), val(d), file(fasta_decoy)

    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_idx.idXML")

    script:
    """
    PeptideIndexer -in ${id_file} \\
    -out ${Sample}_${Condition}_${id}_idx.idXML \\
    -threads ${task.cpus} \\
    -fasta ${fasta_decoy} \\
    -decoy_string DECOY \\
    -enzyme:specificity none
    """
}

