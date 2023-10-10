process OPENMS_PEPTIDEINDEXER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::openms=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.0.0--h8964181_1' :
        'biocontainers/openms:3.0.0--h8964181_1' }"

    input:
        tuple val(meta), path(idxml), path(fasta)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}_${meta.sample}_${meta.condition}_idx"

        """
        PeptideIndexer -in $idxml \\
                -out ${prefix}.idXML \\
                -threads $task.cpus \\
                -fasta $fasta \\
                -decoy_string DECOY \\
                -enzyme:specificity none

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
