process SUMMARIZE_RESULTS {

    conda "bioconda::pyopenms=3.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.3.0--py313h9b5bd11_0' :
        'biocontainers/pyopenms:3.3.0--py313h9b5bd11_0' }"

    input:
    tuple val(meta), path(file)

    output:
    path '*_histogram_mz.csv'           , emit: hist_mz, optional: true
    path '*_histogram_rt.csv'           , emit: hist_rt, optional: true
    path '*_histogram_scores.csv'       , emit: hist_scores, optional: true
    path '*_histogram_xcorr_scores.csv' , emit: hist_xcorr, optional: true
    path '*_peptide_length.csv'         , emit: lengths, optional: true
    path '*_general_stats.csv'          , emit: stats
    path '*.tsv'                        , emit: final_tsv
    path 'versions.yml'                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def quantify = params.quantify ? '--quantify' : ''

    """
    summarize_results.py \\
        --input $file \\
        --out_prefix $prefix \\
        $quantify \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_histogram_mz.csv
    touch ${prefix}_histogram_rt.csv
    touch ${prefix}_histogram_scores.csv
    touch ${prefix}_histogram_xcorr_scores.csv
    touch ${prefix}_peptide_length.csv
    touch ${prefix}_general_stats.csv
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
