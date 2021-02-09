// Import generic module functions
include {  saveFiles } from './functions'

params.options = [:]

process PEAK_PICKING {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), file(mzml_unpicked)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${mzml_unpicked.baseName}.mzML")

    when:
        params.run_centroidisation

    script:
    """
        PeakPickerHiRes -in ${mzml_unpicked} \\
            -out ${mzml_unpicked.baseName}.mzML \\
            -algorithm:ms_levels ${params.pick_ms_levels}
    """
}