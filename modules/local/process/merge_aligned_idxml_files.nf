// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process MERGE_ALIGNED_IDMXL_FILES {

    input:
    tuple val(id), val(Sample), val(Condition), file(ids_aligned)

    output:
    tuple val("$id"), val("$Sample"), val(Condition), file("${Sample}_all_ids_merged.idXML")

    script:
    """
    IDMerger -in $ids_aligned \\
    -out ${Sample}_all_ids_merged.idXML \\
    -threads ${task.cpus}  \\
    -annotate_file_origin  \\
    -merge_proteins_add_PSMs
    """
}