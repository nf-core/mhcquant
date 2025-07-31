process EPICORE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::epicore=0.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epicore:0.1.5--pyhdfd78af_0' :
        'biocontainers/epicore:0.1.5--pyhdfd78af_0' }"

    input:
        path(fasta)
        tuple val(meta), path(result_tsv), path(general_stats)

    output:
        path "*_general_stats.csv",                 emit: stats
        path "${meta.id}.csv",                      emit: final_epicore_tsv
        path "epicore_length_distribution.html",    emit: length_dist
        path "epicore_intensity_histogram.html",    emit: intensity_hist
        path "versions.yml",                        emit: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """#!/bin/bash

        epicore --reference_proteome $fasta --out_dir . generate-epicore-csv $args --evidence_file $result_tsv --html

        mv pep_cores_mapping.csv ${prefix}.csv
        mv length_distributions.html epicore_length_distribution.html
        mv epitope_intensity_hist.html epicore_intensity_histogram.html

        wc -l < epitopes.csv | awk '{print \$1 - 1}' > epicores.txt

        awk 'NR==1 {print \$0 ",# Epicores"; next} NR==2 {getline extra < "epicores.txt"; print \$0 "," extra}' $general_stats > _modified_$general_stats
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            epicore: \$(epicore --version | grep 'epicore' | cut -d ' ' -f3)
        """

    stub:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_general_stats.csv
        touch ${meta.id}_epicore.csv
        touch epicore_length_distribution.html
        touch epicore_intensity_hist.html
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            epicore: \$(epicore --version | grep 'epicore' | cut -d ' ' -f3)
        END_VERSIONS
        """
}
