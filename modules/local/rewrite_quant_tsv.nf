process REWRITE_QUANT_TSV {
    tag "$meta.id"

    conda "bioconda::r-base=4.2.1 bioconda::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c4e658267db06cfc91d5b54f9eb8448b39f0f9d6:12f178b2697c0fd75af5738b089bcdab5cb5fd31-0' :
        'biocontainers/mulled-v2-c4e658267db06cfc91d5b54f9eb8448b39f0f9d6:12f178b2697c0fd75af5738b089bcdab5cb5fd31-0' }"

    input:
        tuple val(meta), path(file)

    output:
        tuple val(meta), path("*.tsv") , emit: mzml
        path "versions.yml"            , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${file.baseName}_final"

        """
        convert_tsv_files.R --input $file --meta_file ${params.input} --outname ${prefix}.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$(echo \$(R --version 2>&1) | awk -F' '  '{print \$3}')
        END_VERSIONS
        """
}
