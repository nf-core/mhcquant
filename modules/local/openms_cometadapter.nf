// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process OPENMS_COMETADAPTER {
    //tag "$meta.id"
    tag "${Sample}"
    label 'process_medium_long'
    
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        // tuple val(meta), path(mzml), path(decoy)
        tuple val(Sample), val(id), val(Condition), path(mzml_file), val(d), path(fasta_decoy)
    
    output:
        // tuple val(meta), path("*.idXML"), emit: idxml
        tuple val("$id"), val("$Sample"), val("$Condition"), path("${Sample}_${Condition}_${id}.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        // def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        def prefix   = options.suffix ? "${Sample}_${options.suffix}" : "${Sample}_${Condition}_${id}"
        
        """
            CometAdapter -in ${mzml_file} \\
                -out ${prefix}.idXML \\
                -database ${fasta_decoy} \\
                -threads $task.cpus $options.args
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """

        // """
        //     CometAdapter  -in ${mzml} -out ${prefix}.idXML -database ${decoy} -threads $task.cpus $options.args
        //     echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        // """

}