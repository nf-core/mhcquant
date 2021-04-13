// Import generic module functions
include {  initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process OPENMS_COMETADAPTER {
    tag "${Sample}"
    label 'process_medium_long'
    
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(Sample), val(id), val(Condition), file(mzml_file), val(d), file(fasta_decoy)
    
    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        
        """
            CometAdapter -in ${mzml_file} \\
            -out ${Sample}_${Condition}_${id}.idXML \\
            -database ${fasta_decoy} \\
            -threads $task.cpus $options.args
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """

}