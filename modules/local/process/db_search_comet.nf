// Import generic module functions
include {  initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

// TODO: Check if the params need to be handed down to the process too

process DB_SEARCH_COMET {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }
    
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(Sample), val(id), val(Condition), file(mzml_file), val(d), file(fasta_decoy)
        val a_ions
        val c_ions
        val x_ions
        val z_ions
        val NL_ions
        val rm_precursor
    
    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}.idXML")

    script:
        """
        CometAdapter  -in ${mzml_file} \\
            -out ${Sample}_${Condition}_${id}.idXML \\
            -threads ${task.cpus} \\
            -database ${fasta_decoy} \\
            -precursor_mass_tolerance ${params.precursor_mass_tolerance} \\
            -fragment_bin_tolerance ${params.fragment_mass_tolerance} \\
            -fragment_bin_offset ${params.fragment_bin_offset} \\
            -num_hits ${params.num_hits} \\
            -digest_mass_range ${params.digest_mass_range} \\
            -max_variable_mods_in_peptide ${params.number_mods} \\
            -allowed_missed_cleavages 0 \\
            -precursor_charge ${params.prec_charge} \\
            -activation_method ${params.activation_method} \\
            -variable_modifications ${params.variable_mods.tokenize(',').collect { "'${it}'" }.join(" ") } \\
            -fixed_modifications ${params.fixed_mods.tokenize(',').collect { "'${it}'"}.join(" ")} \\
            -enzyme '${params.enzyme}' \\
            -spectrum_batch_size ${params.spectrum_batch_size} \\
            $a_ions \\
            $c_ions \\
            $x_ions \\
            $z_ions \\
            $NL_ions \\
            $rm_precursor \\
        """
}