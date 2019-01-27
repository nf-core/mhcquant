#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mhcquant
========================================================================================
 nf-core/mhcquant Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/mhcquant
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/mhcquant : v${workflow.manifest.version}
    =======================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/mhcquant --mzmls '*.mzML' --fasta '*.fasta' --vcf '*.vcf' --alleles 'alleles.tsv' --include_proteins_from_vcf --run_prediction -profile standard,docker

    Mandatory arguments:
      --mzmls                           Path to input data (must be surrounded with quotes)
      --fasta                           Path to Fasta reference
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    Mass Spectrometry Search:
      --peptide_min_length              Minimum peptide length for filtering
      --peptide_max_length              Maximum peptide length for filtering
      --precursor_mass_tolerance        Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance         Mass tolerance of fragment mass bin (ppm)
      --fragment_bin_offset             Offset of fragment mass bin (Comet specific parameter)
      --fdr_threshold                   Threshold for FDR filtering
      --fdr_level                       Level of FDR calculation ('peptide-level-fdrs', 'psm-level-fdrs', 'protein-level-fdrs')
      --digest_mass_range               Mass range of peptides considered for matching
      --activation_method               Fragmentation method ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')
      --enzyme                          Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see OpenMS enzymes)
      --number_mods                     Maximum number of modifications of PSMs
      --fixed_mods                      Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods                   Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --num_hits                        Number of reported hits
      --run_centroidisation             Specify whether mzml data is peak picked or not (true, false)
      --pick_ms_levels                  The ms level used for peak picking (eg. 1, 2)
      --prec_charge                     Precursor charge (eg. "2:3")
      --spectrum_batch_size             Size of Spectrum batch for Comet processing (Decrease/Increase depending on Memory Availability)

    Binding Predictions:
      --run_prediction                  Whether a affinity prediction using MHCFlurry should be run on the results - check if alleles are supported (true, false)
      --refine_fdr_on_predicted_subset  Whether affinity predictions using MHCFlurry should be used to subset PSMs and refine the FDR (true, false)
      --subset_affinity_threshold       Predicted affinity threshold (nM) which will be applied to subset PSMs in FDR refinement. (eg. 500)
      --alleles                         Path to file including allele information

    Variants:
      --include_proteins_from_vcf       Whether to use a provided vcf file to generate proteins and include them in the database search (true, false)
      --vcf                             Path to vcf file
      --variant_annotation_style        Specify which software style was used to carry out the variant annotation in the vcf ("SNPEFF","VEP","ANNOVAR")
      --variant_reference               Specify reference genome used for variant annotation ("GRCH37","GRCH38")
      --variant_indel_filter            Remove insertions and deletions from vcf (true, false)
      --variant_frameshift_filter       Remove insertions and deltionns causing frameshifts from vcf (true, false)
      --variant_snp_filter              Remove snps from vcf (true, false)

    Other options:
      --outdir                          The output directory where the results will be saved
      --email                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                        The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                       The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


// Validate inputs
params.mzmls = params.mzmls ?: { log.error "No read data privided. Make sure you have used the '--mzmls' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No read data privided. Make sure you have used the '--fasta' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()



/*
 * Define the default parameters
 */

//MS params
params.peptide_min_length = 8
params.peptide_max_length = 12
params.fragment_mass_tolerance = 0.02
params.precursor_mass_tolerance = 5
params.fragment_bin_offset = 0
params.fdr_threshold = 0.01
params.fdr_level = 'peptide-level-fdrs'
fdr_level = '-' + params.fdr_level
params.number_mods = 3

params.num_hits = 1
params.digest_mass_range = "800:2500"
params.pick_ms_levels = 2
params.run_centroidisation = false

params.prec_charge = '2:3'
params.activation_method = 'ALL'

params.enzyme = 'unspecific cleavage'
params.fixed_mods = ''
params.variable_mods = 'Oxidation (M)'
params.spectrum_batch_size = 500

//prediction params
params.run_prediction = false
params.refine_fdr_on_predicted_subset = false
params.subset_affinity_threshold = 500

//variant params
params.inlude_proteins_from_vcf = false
params.variant_annotation_style = "SNPEFF"
params.variant_reference = "GRCH38"
params.variant_indel_filter = false
if (params.variant_indel_filter) {
variant_indel_filter="-fINDEL"
} else {
variant_indel_filter=""
}
params.variant_frameshift_filter = false
if (params.variant_frameshift_filter) {
variant_frameshift_filter="-fFS"
} else {
variant_frameshift_filter=""
}
params.variant_snp_filter = false
if (params.variant_snp_filter) {
variant_snp_filter="-fSNP"
} else {
variant_snp_filter=""
}
/*
 * SET UP CONFIGURATION VARIABLES
 */


// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

output_docs = file("$baseDir/docs/output.md")


// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}


/*
 * Create a channel for input mzml files
 */
if( params.run_centroidisation) {
    Channel
        .fromPath( params.mzmls )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.mzmls}\nNB: Path needs to be enclosed in quotes!" }
        .set { input_mzmls_unpicked }

    input_mzmls = Channel.empty()
    input_mzmls_align = Channel.empty()

} else {
    Channel
        .fromPath( params.mzmls )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.mzmls}\nNB: Path needs to be enclosed in quotes!" }
        .into { input_mzmls; input_mzmls_align }

    input_mzmls_unpicked = Channel.empty()
    input_mzmls_align_unpicked = Channel.empty()
}


/*
 * Create a channel for input fasta file
 */
if( params.include_proteins_from_vcf) {
    Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta_vcf}

    input_fasta = Channel.empty()

} else {
    Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta}

    input_fasta_vcf = Channel.empty()
}


/*
 * Create a channel for input alleles file
 */
if( params.run_prediction){
    Channel
        .fromPath( params.alleles )
        .ifEmpty { exit 1, "params.alleles was empty - no input file supplied" }
        .into { input_alleles; input_alleles_refine}
} else {

    input_alleles = Channel.empty()
    input_alleles_refine = Channel.empty()
}




/*
 * Create a channel for input vcf file
 */
if( params.include_proteins_from_vcf){
    Channel
        .fromPath( params.vcf )
        .ifEmpty { exit 1, "params.vcf was empty - no input file supplied" }
        .set { input_vcf}
} else {

    input_vcf = Channel.empty()
}


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/mhcquant v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/mhcquant'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['mzMLs']        = params.mzmls
summary['Fasta Ref']    = params.fasta
summary['Predictions']  = params.run_prediction
summary['SubsetFDR']    = params.refine_fdr_on_predicted_subset
summary['Alleles']      = params.alleles
summary['Variants']     = params.vcf
summary['Centroidisation'] = params.run_centroidisation
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-mhcquant-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/mhcquant Workflow Summary'
    section_href: 'https://github.com/nf-core/mhcquant'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * STEP 0 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * STEP 0.5 - If specified translate variants to proteins and include in reference fasta
 */
process generate_proteins_from_vcf {
    publishDir "${params.outdir}/"

    input:
     file fasta_file_vcf from input_fasta_vcf
     file vcf_file from input_vcf

    output:
     file "${fasta_file_vcf.baseName}_added_vcf.fasta" into appended_fasta

    when:
     params.include_proteins_from_vcf

    script:
     """
     variants2fasta.py -v ${vcf_file} -t ${params.variant_annotation_style} -r ${params.variant_reference} -f ${fasta_file_vcf} -o ${fasta_file_vcf.baseName}_added_vcf.fasta ${variant_indel_filter} ${variant_snp_filter} ${variant_frameshift_filter}
     """
}


/*
 * STEP 1 - generate reversed decoy database
 */
process generate_decoy_database {

    input:
     file fastafile from input_fasta.mix(appended_fasta)

    output:
     file "${fastafile.baseName}_decoy.fasta" into (fastafile_decoy_1, fastafile_decoy_2)
     
    script:
     """
     DecoyDatabase  -in ${fastafile} \\
                    -out ${fastafile.baseName}_decoy.fasta \\
                    -decoy_string DECOY_ \\
                    -decoy_string_position prefix
     """
}


/*
 * STEP 1.5 - Optional: Run Peak Picking as Preprocessing
 */
process peak_picking {

    input:
     file mzml_unpicked from input_mzmls_unpicked

    output:
     file "${mzml_unpicked.baseName}.mzML" into (input_mzmls_picked, input_mzmls_align_picked)

    when:
     params.run_centroidisation

    script:
     """
     PeakPickerHiRes -in ${mzml_unpicked} \\
                     -out ${mzml_unpicked.baseName}.mzML \\
                     -algorithm:ms_levels ${params.pick_ms_levels}
     """
}


/*
 * STEP 2 - run comet database search
 */
process db_search_comet {
 
    input:
     file mzml_file from input_mzmls.mix(input_mzmls_picked)
     file fasta_decoy from fastafile_decoy_1.first()

    output:
     file "${mzml_file.baseName}.idXML" into id_files

    script:
     """
     CometAdapter  -in ${mzml_file} \\
                   -out ${mzml_file.baseName}.idXML \\
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
                   -use_NL_ions true \\
                   -variable_modifications '${params.variable_mods}' \\
                   -fixed_modifications ${params.fixed_mods} \\
                   -enzyme '${params.enzyme}' \\
                   -spectrum_batch_size ${params.spectrum_batch_size}
     """

}


/*
 * STEP 3 - index decoy and target hits
 */
process index_peptides {
 
    input:
     file id_file from id_files
     file fasta_decoy from fastafile_decoy_2.first()

    output:
     file "${id_file.baseName}_idx.idXML" into (id_files_idx, id_files_idx_original)

    script:
     """
     PeptideIndexer -in ${id_file} \\
                    -out ${id_file.baseName}_idx.idXML \\
                    -threads ${task.cpus} \\
                    -fasta ${fasta_decoy} \\
                    -decoy_string DECOY \\
                    -enzyme:specificity none
     """

}


/*
 * STEP 4 - calculate fdr for id based alignment
 */
process calculate_fdr_for_idalignment {
 
    input:
     file id_file_idx from id_files_idx

    output:
     file "${id_file_idx.baseName}_fdr.idXML" into id_files_idx_fdr

    script:
     """
     FalseDiscoveryRate -in ${id_file_idx} \\
                        -out ${id_file_idx.baseName}_fdr.idXML \\
                        -threads ${task.cpus}
     """

}


/*
 * STEP 5 - filter fdr for id based alignment
 */
process filter_fdr_for_idalignment {
 
    input:
     file id_file_idx_fdr from id_files_idx_fdr

    output:
     file "${id_file_idx_fdr.baseName}_filtered.idXML" into id_files_idx_fdr_filtered

    script:
     """
     IDFilter -in ${id_file_idx_fdr} \\
              -out ${id_file_idx_fdr.baseName}_filtered.idXML \\
              -threads ${task.cpus} \\
              -score:pep 0.05  \\
              -remove_decoys
     """

}


/*
 * STEP 6 - compute alignment rt transformation
 */
process align_ids {

    input:
     file id_names from id_files_idx_fdr_filtered.collect{it}

    output:
     file '*.trafoXML' into (id_files_trafo_mzml, id_files_trafo_idxml)

    script:
     def out_names = id_names.collect { it.baseName+'.trafoXML' }.join(' ')
     """
     MapAlignerIdentification -in $id_names \\
                              -trafo_out $out_names
     """

}


/*
 * STEP 7 - align mzML files using trafoXMLs
 */
process align_mzml_files {

    input:
     file id_file_trafo from id_files_trafo_mzml.flatten()
     file mzml_file_align from input_mzmls_align.mix(input_mzmls_align_picked)

    output:
     file "${mzml_file_align.baseName}_aligned.mzML" into mzml_files_aligned

    script:
     """
     MapRTTransformer -in ${mzml_file_align} \\
                      -trafo_in ${id_file_trafo} \\
                      -out ${mzml_file_align.baseName}_aligned.mzML \\
                      -threads ${task.cpus}
     """

}


/*
 * STEP 8 - align unfiltered idXMLfiles using trafoXMLs
 */
process align_idxml_files {

    input:
     file idxml_file_trafo from id_files_trafo_idxml.flatten()
     file idxml_file_align from id_files_idx_original

    output:
     file "${idxml_file_align.baseName}_aligned.idXML" into idxml_files_aligned

    script:
     """
     MapRTTransformer -in ${idxml_file_align} \\
                      -trafo_in ${idxml_file_trafo} \\
                      -out ${idxml_file_align.baseName}_aligned.idXML \\
                      -threads ${task.cpus}
     """

}


/*
 * STEP 9 - merge aligned idXMLfiles
 */
process merge_aligned_idxml_files {

    input:
     file ids_aligned from idxml_files_aligned.collect{it}

    output:
     file "all_ids_merged.idXML" into id_merged
    
    script:
     """
     IDMerger -in $ids_aligned \\
              -out all_ids_merged.idXML \\
              -threads ${task.cpus}  \\
              -annotate_file_origin
     """

}


/*
 * STEP 10 - extract PSM features for Percolator
 */
process extract_psm_features_for_percolator {
    publishDir "${params.outdir}/Intermediate_Results/"
 
    input:
     file id_file_merged from id_merged

    output:
     file "${id_file_merged.baseName}_psm.idXML" into (id_files_merged_psm, id_files_merged_psm_refine, id_files_merged_psm_refine_2)

    script:
     """
     PSMFeatureExtractor -in ${id_file_merged} \\
                         -out ${id_file_merged.baseName}_psm.idXML \\
                         -threads ${task.cpus} 
     """

}


/*
 * STEP 11 - run Percolator
 */
process run_percolator {
    publishDir "${params.outdir}/Intermediate_Results/"
 
    input:
     file id_file_psm from id_files_merged_psm

    output:
     file "${id_file_psm.baseName}_perc.idXML" into id_files_merged_psm_perc

    script:
     """
     PercolatorAdapter -in ${id_file_psm} \\
                       -out ${id_file_psm.baseName}_perc.idXML \\
                       -trainFDR 0.05 \\
                       -testFDR 0.05 \\
                       -threads ${task.cpus} \\
                       -enzyme no_enzyme \\
                       $fdr_level 
     """

}


/*
 * STEP 12 - filter by percolator q-value
 */
process filter_by_q_value {
    publishDir "${params.outdir}/Intermediate_Results/"
 
    input:
     file id_file_perc from id_files_merged_psm_perc

    output:
     file "${id_file_perc.baseName}_filtered.idXML" into id_files_merged_psm_perc_filtered

    when:
     !params.refine_fdr_on_predicted_subset

    script:
     """
     IDFilter -in ${id_file_perc} \\
              -out ${id_file_perc.baseName}_filtered.idXML \\
              -threads ${task.cpus} \\
              -score:pep ${params.fdr_threshold} \\
              -remove_decoys \\
              -length '${params.peptide_min_length}:${params.peptide_max_length}'
     """

}


/*
 * STEP 12.0 - option refine_fdr_on_predicted_subset: filter by percolator q-value
 */
process filter_by_q_value_first {
    publishDir "${params.outdir}/Intermediate_Results/"
    
    input:
     file id_file_perc from id_files_merged_psm_perc
    
    output:
     file "${id_file_perc.baseName}_filtered.idXML" into id_files_merged_psm_perc_filtered_refine

    when:
     params.refine_fdr_on_predicted_subset
    
    script:
     """
     IDFilter -in ${id_file_perc} \\
              -out ${id_file_perc.baseName}_filtered.idXML \\
              -threads ${task.cpus} \\
              -score:pep ${params.fdr_threshold} \\
              -remove_decoys \\
              -length '${params.peptide_min_length}:${params.peptide_max_length}'
     """

}


/*
 * STEP 12.1 - option refine_fdr_on_predicted_subset: export filtered percolator results as mztab
 */
process export_mztab_perc {

    input:
     file percolator_mztab from id_files_merged_psm_perc_filtered_refine

    output:
     file "${percolator_mztab.baseName}.mzTab" into percolator_ids_mztab

    when:
     params.refine_fdr_on_predicted_subset

    script:
     """
     MzTabExporter -in ${percolator_mztab} \\
                   -out ${percolator_mztab.baseName}.mzTab \\
                   -threads ${task.cpus}
     """

}


/*
 * STEP 12.2 - option refine_fdr_on_predicted_subset: export psm results as mztab
 */
process export_mztab_psm {

    input:
     file psm_mztab from id_files_merged_psm_refine

    output:
     file "${psm_mztab.baseName}.mzTab" into psm_ids_mztab

    when:
     params.refine_fdr_on_predicted_subset

    script:
     """
     MzTabExporter -in ${psm_mztab} \\
                   -out ${psm_mztab.baseName}.mzTab \\
                   -threads ${task.cpus}
     """

}


/*
 * STEP 12.3 - option refine_fdr_on_predicted_subset: predict psm results using mhcflurry to shrink search space
 */
process predict_psms {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
     file perc_mztab_file from percolator_ids_mztab
     file psm_mztab_file from psm_ids_mztab
     file allotypes_refine from input_alleles_refine

    output:
     file "peptide_filter.idXML" into peptide_filter

    when:
     params.refine_fdr_on_predicted_subset

    script:
     """
     mhcflurry-downloads fetch models_class1
     mhcflurry_predict_mztab_for_filtering.py ${params.subset_affinity_threshold} ${allotypes_refine} ${perc_mztab_file} ${psm_mztab_file} peptide_filter.idXML
     """
}


/*
 * STEP 12.4 - option refine_fdr_on_predicted_subset: filter psm results by shrinked search space
 */
process filter_psms_by_predictions {
    publishDir "${params.outdir}/Intermediate_Results/"
    
    input:
     file id_file_psm_filtered from id_files_merged_psm_refine_2
     file peptide_filter_file from peptide_filter

    output:
     file "${id_file_psm_filtered.baseName}_pred_filtered.idXML" into id_files_merged_psm_pred_filtered

    when:
     params.refine_fdr_on_predicted_subset    

    script:
     """
     IDFilter -in ${id_file_psm_filtered} \\
              -out ${id_file_psm_filtered.baseName}_pred_filtered.idXML \\
              -whitelist:ignore_modifications \\
              -whitelist:peptides ${peptide_filter_file}\\
              -threads ${task.cpus} \\
     """

}


/*
 * STEP 12.5 - option refine_fdr_on_predicted_subset: recompute percolator fdr on shrinked search space
 */
process run_percolator_on_predicted_subset {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
     file id_file_psm_subset from id_files_merged_psm_pred_filtered

    output:
     file "${id_file_psm_subset.baseName}_perc.idXML" into id_files_merged_psm_pred_perc

    when:
     params.refine_fdr_on_predicted_subset

    script:
     """
     PercolatorAdapter -in ${id_file_psm_subset} \\
                       -out ${id_file_psm_subset.baseName}_perc.idXML \\
                       -trainFDR 0.05 \\
                       -testFDR 0.05 \\
                       -threads ${task.cpus} \\
                       -enzyme no_enzyme \\
                       $fdr_level
     """

}


/*
 * STEP 12.6 - option refine_fdr_on_predicted_subset: filter results by refined fdr
 */
process filter_refined_q_value {
    publishDir "${params.outdir}/Intermediate_Results/"
     
    input:
     file id_file_perc_pred from id_files_merged_psm_pred_perc
     
    output:
     file "${id_file_perc_pred.baseName}_filtered.idXML" into id_files_merged_psm_pred_perc_filtered

    when:
     params.refine_fdr_on_predicted_subset     

    script:
     """      
     IDFilter -in ${id_file_perc_pred} \\
              -out ${id_file_perc_pred.baseName}_filtered.idXML \\
              -threads ${task.cpus} \\
              -score:pep ${params.fdr_threshold} \\
              -remove_decoys \\
              -length '${params.peptide_min_length}:${params.peptide_max_length}'
     """

}


/*
 * STEP 13 - quantify identifications using targeted feature extraction
 */
process quantify_identifications_targeted {
    publishDir "${params.outdir}/Intermediate_Results/"
 
    input:
     file id_file_quant from id_files_merged_psm_perc_filtered.mix(id_files_merged_psm_pred_perc_filtered).first()
     file mzml_quant from mzml_files_aligned

    output:
     file "${mzml_quant.baseName}.featureXML" into (feature_files, feature_files_2)

    script:
     """
     FeatureFinderIdentification -in ${mzml_quant} \\
                                 -id ${id_file_quant} \\
                                 -out ${mzml_quant.baseName}.featureXML \\
                                 -threads ${task.cpus}
     """

}


/*
 * STEP 14 - link extracted features
 */
process link_extracted_features {

    input:
     file feautres from feature_files.collect{it}

    output:
     file "all_features_merged.consensusXML" into consensus_file
    
    script:
     """
     FeatureLinkerUnlabeledKD -in $feautres \\
                              -out 'all_features_merged.consensusXML' \\
                              -threads ${task.cpus}
     """

}


/*
 * STEP 15 - resolve conflicting ids matching to the same feature
 */
process resolve_conflicts {
 
    input:
     file consensus from consensus_file

    output:
     file "${consensus.baseName}_resolved.consensusXML" into (consensus_file_resolved, consensus_file_resolved_2)

    script:
     """
     IDConflictResolver -in ${consensus} \\
                        -out ${consensus.baseName}_resolved.consensusXML \\
                        -threads ${task.cpus}
     """

}


/*
 * STEP 16 - export all information as text to csv
 */
process export_text {
    publishDir "${params.outdir}/"
 
    input:
     file consensus_resolved from consensus_file_resolved

    output:
     file "${consensus_resolved.baseName}.csv" into consensus_text

    script:
     """
     TextExporter -in ${consensus_resolved} \\
                  -out ${consensus_resolved.baseName}.csv \\
                  -threads ${task.cpus} \\
                  -id:add_hit_metavalues 0 \\
                  -id:add_metavalues 0 \\
                  -id:peptides_only
     """

}


/*
 * STEP 17 - export all information as mzTab
 */
process export_mztab {
    publishDir "${params.outdir}/"

    input:
     file feature_file_2 from consensus_file_resolved_2

    output:
     file "${feature_file_2.baseName}.mzTab" into features_mztab

    script:
     """
     MzTabExporter -in ${feature_file_2} \\
                   -out ${feature_file_2.baseName}.mzTab \\
                   -threads ${task.cpus}
     """

}


/*
 * STEP 18 - If specified predict peptides using MHCFlurry
 */
process predict_peptides {
    publishDir "${params.outdir}/"

    input:
     file mztab_file from features_mztab
     file allotypes from input_alleles

    output:
     file "*predicted_peptides.csv" into predicted_peptides

    when:
     params.run_prediction

    script:
     """
     mhcflurry-downloads fetch models_class1
     mhcflurry_predict_mztab.py ${allotypes} ${mztab_file} predicted_peptides.csv
     """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/mhcquant] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/mhcquant] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/mhcquant] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/mhcquant] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_ir = new File( "${params.outdir}/Intermediate_Results/" )
    if( !output_ir.exists() ) {
      output_ir.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/mhcquant] Pipeline Complete"

}
