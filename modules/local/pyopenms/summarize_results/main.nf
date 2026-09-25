process SUMMARIZE_RESULTS {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.4.1--py312h6b06db6_2' :
        'biocontainers/pyopenms:3.4.1--py312h6b06db6_2' }"

    input:
    tuple val(meta), path(file), path(frag_mass_errs), path(trafoxmls)

    output:
    path '*_histogram_mz.csv'                                   , emit: hist_mz, optional: true
    path '*_histogram_rt.csv'                                   , emit: hist_rt, optional: true
    path '*_histogram_scores.csv'                               , emit: hist_scores, optional: true
    path '*_xcorr_scores.json'                                  , emit: xcorr, optional: true
    path '*_peptide_length.csv'                                 , emit: lengths, optional: true
    path '*_peptide_intensity.json'                             , emit: intensities, optional: true
    path '*_histogram_im.csv'                                   , emit: hist_im, optional: true
    path '*_deeplc_rt_diff.json'                                , emit: rt_calibration, optional: true
    path '*_aligned_residuals.json'                             , emit: aligned_residuals, optional: true
    path '*_frag_mass_err.json'                                 , emit: frag_mass_err, optional: true
    tuple val(meta), path('*.tsv'), path('*_general_stats.csv') , emit: epicore_input
    tuple val("${task.process}"), val('pyopenms'), eval("pip show pyopenms | grep Version | sed 's/Version: //'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def quantify = params.quantify ? '--quantify' : ''
    def frag_mass_err = frag_mass_errs ? "--frag_mass_err ${frag_mass_errs.join(' ')}" : ''
    def trafo = trafoxmls ? "--trafoxml ${trafoxmls.join(' ')}" : ''

    """
    summarize_results.py \\
        --input $file \\
        --out_prefix $prefix \\
        $quantify \\
        $frag_mass_err \\
        $trafo \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_histogram_mz.csv
    touch ${prefix}_histogram_rt.csv
    touch ${prefix}_histogram_scores.csv
    touch ${prefix}_xcorr_scores.json
    touch ${prefix}_peptide_length.csv
    touch ${prefix}_peptide_intensity.json
    touch ${prefix}_histogram_im.csv
    touch ${prefix}_deeplc_rt_diff.json
    touch ${prefix}_aligned_residuals.json
    touch ${prefix}_frag_mass_err.json
    touch ${prefix}_general_stats.csv
    touch ${prefix}.tsv
    """
}
