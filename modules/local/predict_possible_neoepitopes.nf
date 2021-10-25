// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSIONFRED2 = '2.0.6'
def VERSIONMHCNUGGETS = '2.3.2'

process PREDICT_POSSIBLE_NEOEPITOPES {
    tag "$meta"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', publish_id:'') }

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(meta), val(alleles), path(vcf)

    output:
        tuple val(meta), path("${prefix}.csv"), emit: csv
        tuple val(meta), path("${prefix}.txt"), emit: txt
        path "versions.yml", emit: versions

    script:
        def prefix = options.suffix ? "${meta}_${options.suffix}" : "${meta}_vcf_neoepitopes_class1"

        """
            vcf_neoepitope_predictor.py -t ${params.variant_annotation_style} -r ${params.variant_reference} -a '${alleles}' -minl ${params.peptide_min_length} -maxl ${params.peptide_max_length} -v ${vcf} -o ${prefix}.csv

            cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                mhcflurry: \$(mhcflurry-predict --version | sed 's/^mhcflurry //; s/ .*\$//' )
                mhcnuggets: \$(echo $VERSIONMHCNUGGETS)
                FRED2: \$(echo $VERSIONFRED2)
            END_VERSIONS
        """
}
