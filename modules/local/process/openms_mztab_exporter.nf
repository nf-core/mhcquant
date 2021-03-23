// Import generic module functions
include { initOptions; saveFiles } from './functions'

def options    = initOptions(params.options)

//TODO: combine in a subflow --> when needs to be removed
process OPENMS_MZTAB_EXPORTER {
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
        tuple val(id), val(Sample), val(Condition), file(mztab)

    output:
        tuple val("$id"), val("$Sample"), file("*.mzTab"), emit: mztab   
        path  "*.version.txt", emit: version

    script:
    def prefix = options.suffix ? "${Sample}.${options.suffix}" : "${Sample}"

    """
        MzTabExporter -in ${mztab} \\
            -out ${Sample}_${prefix}.mzTab \\
            -threads ${task.cpus}

        FileInfo --help &> openms.version.txt
    """
}