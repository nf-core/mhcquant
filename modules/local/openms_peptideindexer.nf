// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process OPENMS_PEPTIDEINDEXER {

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(Sample), val(id), val(Condition), path(id_file), val(d), path(fasta_decoy)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), path("${Sample}_${Condition}_${id}_idx.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)

        """
            PeptideIndexer -in ${id_file} \\
                    -out ${Sample}_${Condition}_${id}_idx.idXML \\
                    -threads ${task.cpus} \\
                    -fasta ${fasta_decoy} \\
                    -decoy_string DECOY \\
                    -enzyme:specificity none
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}

