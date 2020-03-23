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
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/mhcquant --mzmls '*.mzML' --fasta '*.fasta' --vcf '*.vcf' --class_1_alleles 'alleles.tsv' --include_proteins_from_vcf --predict_class_1 --refine_fdr_on_predicted_subset -profile standard,docker

    Mandatory arguments:
      --mzmls [file]                            Path to input data (must be surrounded with quotes)	
      --raw_input [bool]                        Specify whether raw files should be used as input
      --raw_files [file]                        Path to input raw data (must be surrounded with quotes)
      --fasta [file]                            Path to Fasta reference
      -profile [str]                            Configuration profile to use. Can use multiple (comma separated)
                                                Available: docker, singularity, test, awsbatch and more
    Mass Spectrometry Search:
      --peptide_min_length [int]                Minimum peptide length for filtering
      --peptide_max_length [int]                Maximum peptide length for filtering
      --precursor_mass_tolerance [int           Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance [int]           Mass tolerance of fragment mass bin (ppm)
      --fragment_bin_offset [int]               Offset of fragment mass bin (Comet specific parameter)
      --use_x_ions [bool]                       Use x ions for spectral matching in addition
      --use_z_ions [bool]                       Use z ions for spectral matching in addition
      --use_a_ions [bool]                       Use a ions for spectral matching in addition
      --use_c_ions [bool]                       Use c ions for spectral matching in addition
      --fdr_threshold [int]                     Threshold for FDR filtering
      --fdr_level [str]                         Level of FDR calculation ('peptide-level-fdrs', 'psm-level-fdrs', 'protein-level-fdrs')
      --digest_mass_range [int]                 Mass range of peptides considered for matching
      --activation_method [str]                 Fragmentation method ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')
      --enzyme [str]                            Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see OpenMS enzymes)
      --number_mods [int]                       Maximum number of modifications of PSMs
      --fixed_mods [str]                        Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods [str]                     Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --num_hits [int]                          Number of reported hits
      --run_centroidisation [bool]              Specify whether mzml data is peak picked or not (true, false)
      --pick_ms_levels [int]                    The ms level used for peak picking (eg. 1, 2)
      --prec_charge [str]                       Precursor charge (eg. "2:3")
      --max_rt_alignment_shift [int]            Maximal retention time shift (sec) resulting from linear alignment      
      --spectrum_batch_size [int]               Size of Spectrum batch for Comet processing (Decrease/Increase depending on Memory Availability)
      --description_correct_features [int]      Description of correct features for Percolator (0, 1, 2, 4, 8, see Percolator retention time and calibration) 
      --klammer [bool]                          Retention time features are calculated as in Klammer et al. instead of with Elude.
      --predict_RT [bool]                       Retention time prediction for identified peptides
      --skip_decoy_generation [bool]            Use a fasta database that already includes decoy sequences
      --quantification_fdr [bool]               Assess and assign ids matched between runs with an additional quantification FDR
      --quantification_min_prob  [int]          Specify a minimum probability cut off for quantification

    Binding Predictions:	
      --predict_class_1 [bool]                  Whether a class 1 affinity prediction using MHCFlurry should be run on the results - check if alleles are supported (true, false)	
      --predict_class_2 [bool]                  Whether a class 2 affinity prediction using MHCNuggets should be run on the results - check if alleles are supported (true, false) 	
      --refine_fdr_on_predicted_subset[bool]    Whether affinity predictions using MHCFlurry should be used to subset PSMs and refine the FDR (true, false)	
      --subset_affinity_threshold [int]         Predicted affinity threshold (nM) which will be applied to subset PSMs in FDR refinement. (eg. 500)	
      --class_1_alleles [file]                  Path to file including class 1 allele information	
      --class_2_alleles [file]                  Path to file including class 2 allele information	

    Variants:	
      --include_proteins_from_vcf [bool]        Whether to use a provided vcf file to generate proteins and include them in the database search (true, false)	
      --vcf [file]                              Path to vcf file	
      --variant_annotation_style [str]          Specify which software style was used to carry out the variant annotation in the vcf ("SNPEFF","VEP","ANNOVAR")	
      --variant_reference [str]                 Specify reference genome used for variant annotation ("GRCH37","GRCH38")	
      --variant_indel_filter [bool]             Remove insertions and deletions from vcf (true, false)	
      --variant_frameshift_filter [bool]        Remove insertions and deltionns causing frameshifts from vcf (true, false)	
      --variant_snp_filter [bool]               Remove snps from vcf (true, false)

    Other options:
      --outdir [file]                           The output directory where the results will be saved
      --email [email]                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]                   Same as --email, except only send mail if the workflow is not successful
      -name [str]                               Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                          The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]                         The AWS Region for your AWS Batch job to run on
      --awscli [str]                            Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if(!params.raw_input){
   params.mzmls = params.mzmls ?: { log.error "No spectra data privided. Make sure you have used the '--mzmls' option."; exit 1 }()
} else {
   params.raw_files = params.raw_files ?: { log.error "No spectra data privided. Make sure you have used the '--raw_files' option."; exit 1 }()
}
params.fasta = params.fasta ?: { log.error "No database fasta privided. Make sure you have used the '--fasta' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()

/*
 * SET UP CONFIGURATION VARIABLES
 */

//MS params
params.peptide_min_length = 8
params.peptide_max_length = 12
params.fragment_mass_tolerance = 0.02
params.precursor_mass_tolerance = 5
params.use_x_ions = false
x_ions = params.use_x_ions ? '-use_X_ions true' : ''
params.use_z_ions = false
z_ions = params.use_z_ions ? '-use_Z_ions true' : ''
params.use_a_ions = false
a_ions = params.use_a_ions ? '-use_A_ions true' : ''
params.use_c_ions = false
c_ions = params.use_c_ions ? '-use_C_ions true' : ''
params.fragment_bin_offset = 0
params.fdr_threshold = 0.01
params.fdr_level = 'peptide-level-fdrs'
fdr_level = (params.fdr_level == 'psm-level-fdrs') ? '' : '-'+params.fdr_level
params.description_correct_features = 0
params.klammer = false
params.predict_RT = false
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

params.skip_decoy_generation = false
if (params.skip_decoy_generation) {
log.warn "Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement"
log.warn "Decoys have to be named with DECOY_ as prefix in your fasta database"
}

params.quantification_fdr = false
params.quantification_min_prob = 0
if (params.quantification_fdr) {
   log.warn "Quantification FDR enabled"
}

//prediction params
params.predict_class_1 = false
params.predict_class_2 = false
params.refine_fdr_on_predicted_subset = false
if (params.skip_decoy_generation) {
log.warn "Be aware: subset FDR refinement only considers MHC class I alleles supported by mhcflurry"
}
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

if (params.raw_input){

  input_mzmls = Channel.empty()
  input_mzmls_align = Channel.empty()
  input_mzmls_unpicked = Channel.empty()
  input_mzmls_align_unpicked = Channel.empty()

  Channel
        .fromPath( params.raw_files )
        .ifEmpty { exit 1, "Cannot find any raw files matching: ${params.raw_files}\nNB: Path needs to be enclosed in quotes!" }
        .set { input_raws }


} else {

  input_raws = Channel.empty()

  if (params.run_centroidisation) {
    Channel
        .fromPath( params.mzmls )
        .ifEmpty { exit 1, "Cannot find any mzmls matching: ${params.mzmls}\nNB: Path needs to be enclosed in quotes!" }
        .set { input_mzmls_unpicked }

    input_mzmls = Channel.empty()
    input_mzmls_align = Channel.empty()

  } else {
    Channel
        .fromPath( params.mzmls )
        .ifEmpty { exit 1, "Cannot find any mzmls matching: ${params.mzmls}\nNB: Path needs to be enclosed in quotes!" }
        .into { input_mzmls; input_mzmls_align }

    input_mzmls_unpicked = Channel.empty()
    input_mzmls_align_unpicked = Channel.empty()
  }

}

/*
 * Create a channel for input fasta file
 */
if( params.include_proteins_from_vcf) {
    Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta_vcf }

    input_fasta = Channel.empty()
    input_fasta_1 = Channel.empty()
    input_fasta_2 = Channel.empty()

} else if( params.skip_decoy_generation) {
    Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .into { input_fasta; input_fasta_1; input_fasta_2 }

    input_fasta_vcf = Channel.empty()

} else {
    Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta }

    input_fasta_vcf = Channel.empty()
    input_fasta_1 = Channel.empty()
    input_fasta_2 = Channel.empty()

}


/*
 * Create a channel for class 1 alleles file
 */
if( params.predict_class_1){
    Channel
        .fromPath( params.class_1_alleles )
        .ifEmpty { exit 1, "params.alleles was empty - no input file supplied" }
        .into { peptides_class_1_alleles; peptides_class_1_alleles_refine; neoepitopes_class_1_alleles; neoepitopes_class_1_alleles_prediction}

} else {

    peptides_class_1_alleles = Channel.empty()
    peptides_class_1_alleles_refine = Channel.empty()
    neoepitopes_class_1_alleles = Channel.empty()
    neoepitopes_class_1_alleles_prediction = Channel.empty()
}

/*
 * Create a channel for class  2 allele files
*/
if( params.predict_class_2){
    Channel
        .fromPath( params.class_2_alleles )
        .ifEmpty { exit 1, "params.class_2_alleles was empty - no input file supplied" }
        .into { nepepitopes_class_2_alleles; peptides_class_2_alleles; peptides_class_2_alleles_II }

} else {
    nepepitopes_class_2_alleles = Channel.empty()
    peptides_class_2_alleles = Channel.empty()
    peptides_class_2_alleles_II = Channel.empty()
}

/*
 * Create a channel for input vcf file
 */
if( params.include_proteins_from_vcf){
    Channel
        .fromPath( params.vcf )
        .ifEmpty { exit 1, "params.vcf was empty - no input file supplied" }
        .into { input_vcf; input_vcf_neoepitope; input_vcf_neoepitope_II}

} else {
    input_vcf = Channel.empty()
    input_vcf_neoepitope = Channel.empty()
    input_vcf_neoepitope_II = Channel.empty()
}

// Stage config files
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']           = custom_runName ?: workflow.runName
summary['Pipeline Name']      = 'nf-core/mhcquant'
summary['Pipeline Version']   = workflow.manifest.version
summary['Run Name']           = custom_runName ?: workflow.runName
summary['mzMLs']              = params.mzmls
summary['Fasta Ref']          = params.fasta
summary['Class 1 Prediction'] = params.predict_class_1
summary['Class 2 Prediction'] = params.predict_class_2
summary['SubsetFDR']          = params.refine_fdr_on_predicted_subset
summary['Quantification FDR'] = params.quantification_fdr
summary['Class 1 Alleles']    = params.predict_class_1
summary['Class 2 Alelles']    = params.predict_class_2
summary['RT Prediction']      = params.predict_RT
summary['Variants']           = params.include_proteins_from_vcf
summary['Centroidisation']    = params.run_centroidisation
summary['Max Memory']         = params.max_memory
summary['Max CPUs']           = params.max_cpus
summary['Max Time']           = params.max_time
summary['Output dir']         = params.outdir
summary['Working dir']        = workflow.workDir
summary['Container Engine']   = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']       = "$HOME"
summary['Current user']       = "$USER"
summary['Current path']       = "$PWD"
summary['Working dir']        = workflow.workDir
summary['Output dir']         = params.outdir
summary['Script dir']         = workflow.projectDir
summary['Config Profile']     = workflow.profile
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']     = params.awsregion
    summary['AWS Queue']      = params.awsqueue
    summary['AWS CLI']        = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-mhcquant-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/mhcquant Workflow Summary'
    section_href: 'https://github.com/nf-core/mhcquant'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    FileInfo --help &> v_openms.txt
    percolator -h &> v_percolator.txt
    comet -p
    mhcflurry-predict --version &> v_mhcflurry.txt
    echo 2.3.1 > v_mhcnuggets.txt
    echo 2.0.3 > v_fred2.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * STEP 0 - Output Description HTML
*/ 
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
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
    
    when:
     !params.skip_decoy_generation
 
    script:
     """
     DecoyDatabase  -in ${fastafile} \\
                    -out ${fastafile.baseName}_decoy.fasta \\
                    -decoy_string DECOY_ \\
                    -decoy_string_position prefix
     """
}


/*
 * STEP 1.1 - Raw file conversion
 */
process raw_file_conversion {

    input:
     file rawfile from input_raws

    output:
     file "${rawfile.baseName}.mzML" into (raws_converted, raws_converted_align)
   
    when:
     params.raw_input
    
    script:
     """
     ThermoRawFileParser.sh -i=${rawfile} -f=2 -b=${rawfile.baseName}.mzML
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

    label 'process_high'
 
    input:
     file mzml_file from raws_converted.mix(input_mzmls.mix(input_mzmls_picked))
     file fasta_decoy from fastafile_decoy_1.mix(input_fasta_1).first()

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
                   -variable_modifications ${params.variable_mods.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                   -fixed_modifications ${params.fixed_mods.tokenize(',').collect { "'${it}'"}.join(" ")} \\
                   -enzyme '${params.enzyme}' \\
                   -spectrum_batch_size ${params.spectrum_batch_size} \\
                   $a_ions \\
                   $c_ions \\
                   $x_ions \\
                   $z_ions \\     
     """

}


/*
 * STEP 3 - index decoy and target hits
 */
process index_peptides {
 
    input:
     file id_file from id_files
     file fasta_decoy from fastafile_decoy_2.mix(input_fasta_2).first()

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
     file "${id_file_idx_fdr.baseName}_filtered.idXML" into (id_files_idx_fdr_filtered, id_files_for_quant_fdr)

    script:
     """
     IDFilter -in ${id_file_idx_fdr} \\
              -out ${id_file_idx_fdr.baseName}_filtered.idXML \\
              -threads ${task.cpus} \\
              -score:pep ${params.fdr_threshold} \\
              -length '${params.peptide_min_length}:${params.peptide_max_length}' \\
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
     file '*.trafoXML' into id_files_trafo

    script:
     def out_names = id_names.collect { it.baseName+'.trafoXML' }.join(' ')
     """
     MapAlignerIdentification -in $id_names \\
                              -trafo_out $out_names \\
                              -model:type linear \\
                              -algorithm:max_rt_shift ${params.max_rt_alignment_shift}
     """

}

input_mzmls_align
 .mix(raws_converted_align)
 .mix(input_mzmls_align_picked)
 .collectFile( sort: { it.baseName } )
 .set{input_mzmls_combined}

id_files_idx_original
 .collectFile( sort: { it.baseName } )
 .set{input_ids_sorted}

id_files_for_quant_fdr
 .collectFile( sort: { it.baseName } )
 .set{input_ids_for_quant_fdr_sorted}

id_files_trafo
 .flatten()
 .collectFile( sort: { it.baseName } )
 .into{trafo_sorted_mzml; trafo_sorted_id}


/*
 * STEP 7 - align mzML files using trafoXMLs
 */
process align_mzml_files {

    input:
     file id_file_trafo from trafo_sorted_mzml
     file mzml_file_align from input_mzmls_combined

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
     file idxml_file_trafo from trafo_sorted_id
     file idxml_file_align from input_ids_sorted

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

    if (params.klammer && params.description_correct_features == 0) {
        log.warn('Klammer was specified, but description of correct features was still 0. Please provide a description of correct features greater than 0.')
        log.warn('Klammer has been turned off!')
    }

    script:
    if (params.description_correct_features > 0 && params.klammer){
    """
    PercolatorAdapter -in ${id_file_psm} \\
                       -out ${id_file_psm.baseName}_perc.idXML \\
                       -trainFDR 0.05 \\
                       -testFDR 0.05 \\
                       -threads ${task.cpus} \\
                       -enzyme no_enzyme \\
                       $fdr_level \\
                       -doc ${params.description_correct_features} \\
                       -klammer
    """
    } else {
    """
    PercolatorAdapter -in ${id_file_psm} \\
                       -out ${id_file_psm.baseName}_perc.idXML \\
                       -trainFDR 0.05 \\
                       -testFDR 0.05 \\
                       -threads ${task.cpus} \\
                       -enzyme no_enzyme \\
                       $fdr_level \\
                       -doc ${params.description_correct_features} \\
    """
    }
     

}


/*
 * STEP 12 - filter by percolator q-value
 */
process filter_by_q_value {
    publishDir "${params.outdir}/Intermediate_Results/"
 
    input:
     file id_file_perc from id_files_merged_psm_perc

    output:
     file "${id_file_perc.baseName}_filtered.idXML" into (id_files_merged_psm_perc_filtered, ids_for_rt_training, ids_for_rt_prediction)

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
    publishDir "${params.outdir}/Intermediate_Results/"

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
    publishDir "${params.outdir}/Intermediate_Results/"

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
    echo true

    input:
     file perc_mztab_file from percolator_ids_mztab
     file psm_mztab_file from psm_ids_mztab
     file allotypes_refine from peptides_class_1_alleles_refine

    output:
     file "peptide_filter.idXML" into peptide_filter

    when:
     params.refine_fdr_on_predicted_subset

    script:
     """
     mhcflurry-downloads --quiet fetch models_class1
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
     file "${id_file_perc_pred.baseName}_filtered.idXML" into (id_files_merged_psm_pred_perc_filtered, ids_for_rt_training_subset, ids_for_rt_prediction_subset)

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
     file id_file_quant_int from input_ids_for_quant_fdr_sorted

    output:
     file "${mzml_quant.baseName}.featureXML" into (feature_files, feature_files_II)

    script:
    if (!params.quantification_fdr){
     """
     FeatureFinderIdentification -in ${mzml_quant} \\
                                 -id ${id_file_quant} \\
                                 -out ${mzml_quant.baseName}.featureXML \\
                                 -threads ${task.cpus}
     """
    } else {
          """
     FeatureFinderIdentification -in ${mzml_quant} \\
                                 -id ${id_file_quant_int} \\
                                 -id_ext ${id_file_quant} \\
                                 -svm:min_prob ${params.quantification_min_prob} \\
                                 -out ${mzml_quant.baseName}.featureXML \\
                                 -threads ${task.cpus}
     """   
    }
}


/*
 * STEP 14 - link extracted features
 */
process link_extracted_features {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
     file features from feature_files.collect{it}
     file features_base from feature_files_II.map { file -> file.baseName + '\n' }.collect{it}

    output:
     file "mhcquant_file_order.txt" into order_file
     file "all_features_merged.consensusXML" into consensus_file
    
    script:
     """
     cat ${features_base} > 'mhcquant_file_order.txt'
     FeatureLinkerUnlabeledKD -in ${features} \\
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
     file "${feature_file_2.baseName}.mzTab" into features_mztab, features_mztab_neoepitopes, features_mztab_neoepitopes_II, mhcnuggets_mztab

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
process predict_peptides_mhcflurry_class_1 {
    publishDir "${params.outdir}/class_1_bindings"
    echo true

    input:
     file mztab_file from features_mztab
     file class_1_alleles from peptides_class_1_alleles

    output:
     file "*predicted_peptides_class_1.csv" into predicted_peptides

    when:
     params.predict_class_1

    script:
     """
     mhcflurry-downloads --quiet fetch models_class1
     mhcflurry_predict_mztab.py ${class_1_alleles} ${mztab_file} predicted_peptides_class_1.csv
     """
}


/*
 * STEP 19 - Preprocess found peptides for MHCNuggets prediction class 2
 */ 
 process preprocess_peptides_mhcnuggets_class_2 {
     
    input:
     file mztab_file from mhcnuggets_mztab

    output:
     file 'preprocessed_mhcnuggets_peptides' into preprocessed_mhcnuggets_peptides
     file 'peptide_to_geneID' into peptide_to_geneID

    when:
     params.predict_class_2

    script:
    """
    preprocess_peptides_mhcnuggets.py --mztab ${mztab_file} --output preprocessed_mhcnuggets_peptides
    """
 }


 /*
 * STEP 20 - Predict found peptides using MHCNuggets class 2
*/  
 process predict_peptides_mhcnuggets_class_2 {

    input:
     file preprocessed_peptides from preprocessed_mhcnuggets_peptides
     file class_2_alleles from peptides_class_2_alleles

    output:
     file '*_predicted_peptides_class_2' into predicted_mhcnuggets_peptides

    when:
     params.predict_class_2

    script:
    """
    mhcnuggets_predict_peptides.py --peptides ${preprocessed_peptides} --alleles ${class_2_alleles} --output _predicted_peptides_class_2
    """
 }


 /*
 * STEP 21 - Postprocess predicted MHCNuggets peptides class 2
 */ 
 process postprocess_peptides_mhcnuggets_class_2 {
    publishDir "${params.outdir}/class_2_bindings"

    input:
     file predicted_peptides from predicted_mhcnuggets_peptides.collect{it}
     file peptide_to_geneID from peptide_to_geneID

    output:
     file '*.csv' into postprocessed_peptides_mhcnuggets

    when:
     params.predict_class_2

    script:
    """
    postprocess_peptides_mhcnuggets.py --input ${predicted_peptides} --peptides_seq_ID ${peptide_to_geneID}
    """
 }
 

/*
 * STEP 22 - Predict all possible neoepitopes from vcf
 */
process predict_possible_neoepitopes {
    publishDir "${params.outdir}/"
    echo true

    label 'process_high'

    input:
     file alleles_file from neoepitopes_class_1_alleles
     file vcf_file from input_vcf_neoepitope

    output:
     file "vcf_neoepitopes.csv" into possible_neoepitopes
     file "vcf_neoepitopes.txt" into possible_neoepitopes_list
 
    when:
     params.include_proteins_from_vcf
     params.predict_class_1

    script:
     """
     vcf_neoepitope_predictor.py -t ${params.variant_annotation_style} -r ${params.variant_reference} -a ${alleles_file} -minl ${params.peptide_min_length} -maxl ${params.peptide_max_length} -v ${vcf_file} -o vcf_neoepitopes.csv
     """
}


/*
 * STEP 22/2 - Predict all possible neoepitopes from vcf
 */
process predict_possible_class_2_neoepitopes {
    publishDir "${params.outdir}/"
    echo true

    label 'process_high'

    input:
     file alleles_file_II from peptides_class_2_alleles_II
     file vcf_file from input_vcf_neoepitope_II

    output:
     file "vcf_neoepitopes.csv" into possible_neoepitopes_II
     file "vcf_neoepitopes.txt" into possible_neoepitopes_list_II

    when:
     params.include_proteins_from_vcf
     !params.predict_class_1
     params.predict_class_2

    script:
     """
     vcf_neoepitope_predictor.py -t ${params.variant_annotation_style} -r ${params.variant_reference} -a ${alleles_file_II} -minl ${params.peptide_min_length} -maxl ${params.peptide_max_length} -v ${vcf_file} -o vcf_neoepitopes.csv
     """
}


/*
 * STEP 23 - Resolve found neoepitopes
 */
process Resolve_found_neoepitopes {
    publishDir "${params.outdir}/"
    echo true

    input:
     file mztab from features_mztab_neoepitopes
     file neoepitopes from possible_neoepitopes

    output:
     file "found_neoepitopes_class_1.csv" into found_neoepitopes
    
    when:
     params.include_proteins_from_vcf
     params.predict_class_1

    script:
     """
     resolve_neoepitopes.py -n ${neoepitopes} -m ${mztab} -f csv -o found_neoepitopes_class_1
     """
}


/*
 * STEP 23/2 - Resolve found neoepitopes
 */
process Resolve_found_class_2_neoepitopes {
    publishDir "${params.outdir}/"
    echo true

    input:
     file mztab from features_mztab_neoepitopes_II
     file neoepitopes from possible_neoepitopes_II

    output:
     file "found_neoepitopes_class_2.csv" into found_neoepitopes_II, mhcnuggets_neo_preprocessing, mhcnuggets_neo_postprocessing

    when:
     params.include_proteins_from_vcf
     params.predict_class_2

    script:
     """
     resolve_neoepitopes.py -n ${neoepitopes} -m ${mztab} -f csv -o found_neoepitopes_class_2
     """
}


/*
 * STEP 24 - Predict class 1 neoepitopes MHCFlurry
 */
process Predict_neoepitopes_mhcflurry_class_1 {
    publishDir "${params.outdir}/class_1_bindings"
    echo true

    input:
     file allotypes from neoepitopes_class_1_alleles_prediction
     file neoepitopes from found_neoepitopes

    output:
     file "*predicted_neoepitopes_class_1.csv" into predicted_neoepitopes
    
    when:
     params.include_proteins_from_vcf
     params.predict_class_1

    script:
     """
     mhcflurry-downloads --quiet fetch models_class1
     mhcflurry_neoepitope_binding_prediction.py ${allotypes} ${neoepitopes} predicted_neoepitopes_class_1.csv
     """
}


/*
 * STEP 25 - Preprocess resolved neoepitopes in a format that MHCNuggets understands
 */
process preprocess_neoepitopes_mhcnuggets_class_2 {

    input:
    file neoepitopes from mhcnuggets_neo_preprocessing

    output:
    file 'mhcnuggets_preprocessed' into preprocessed_mhcnuggets_neoepitopes

    when:
     params.include_proteins_from_vcf
     params.predict_class_2

    script:
    """
    preprocess_neoepitopes_mhcnuggets.py --neoepitopes ${neoepitopes} --output mhcnuggets_preprocessed
    """
}


/*
 * STEP 26 - Predict class 2 MHCNuggets
 */
process predict_neoepitopes_mhcnuggets_class_2 {

    input:
    file preprocessed_neoepitopes from preprocessed_mhcnuggets_neoepitopes
    file cl_2_alleles from nepepitopes_class_2_alleles

    output:
    file '*predicted_neoepitopes_class_2' into predicted_neoepitopes_class_2

    when:
     params.include_proteins_from_vcf
     params.predict_class_2

    script:
    """
    mhcnuggets_predict_peptides.py --peptides ${preprocessed_neoepitopes} --alleles ${cl_2_alleles} --output _predicted_neoepitopes_class_2
    """
}


/*
 * STEP 27 - Class 2 MHCNuggets Postprocessing
*/ 
process postprocess_neoepitopes_mhcnuggets_class_2 {
    publishDir "${params.outdir}/class_2_bindings"

    input:
    file neoepitopes from mhcnuggets_neo_postprocessing
    file predicted_cl_2 from predicted_neoepitopes_class_2.collect{it}

    output:
    file '*.csv' into postprocessed_predicted_neoepitopes_class_2

    when:
     params.include_proteins_from_vcf
     params.predict_class_2

    script:
    """
    postprocess_neoepitopes_mhcnuggets.py --input ${predicted_cl_2} --neoepitopes ${neoepitopes}
    """
}


/*
 * STEP 28 - Train Retention Times Predictor
*/
process train_retention_time_predictor {

    input:
    file id_files_for_rt_training from ids_for_rt_training.mix(ids_for_rt_training_subset)

    output:
    file "${id_files_for_rt_training.baseName}.txt" into (trained_rt_model, trained_rt_model_II)
    file "${id_files_for_rt_training.baseName}_params.paramXML" into (trained_rt_params, trained_rt_params_II) 
    file "${id_files_for_rt_training.baseName}_trainset.txt" into (trained_rt_set,  trained_rt_set_II)

    when:
     params.predict_RT

    script:
    """
    RTModel -in ${id_files_for_rt_training} \\
            -cv:skip_cv \\
            -out "${id_files_for_rt_training.baseName}.txt" \\
            -out_oligo_params "${id_files_for_rt_training.baseName}_params.paramXML" \\
            -out_oligo_trainset "${id_files_for_rt_training.baseName}_trainset.txt"
    """
}


/*
 * STEP 29 - Retention Times Predictor Found Peptides
*/
process predict_retention_times_of_found_peptides {
    publishDir "${params.outdir}/RT_prediction/"

    input:
    file id_files_for_rt_prediction from ids_for_rt_prediction.mix(ids_for_rt_prediction_subset)
    file trained_rt_param_file from trained_rt_params
    file trained_rt_set_file from trained_rt_set
    file trained_rt_model_file from trained_rt_model

    output:
    file "${id_files_for_rt_prediction.baseName}_RTpredicted.csv" into rt_predicted

    when:
     params.predict_RT

    script:
    """
    RTPredict -in_id ${id_files_for_rt_prediction} \\
              -svm_model ${trained_rt_model_file} \\
              -in_oligo_params ${trained_rt_param_file} \\
              -in_oligo_trainset ${trained_rt_set_file} \\
              -out_text:file "${id_files_for_rt_prediction.baseName}_RTpredicted.csv"
    """
}


/*
 * STEP 29 - Retention Times Predictor possible Neoepitopes
*/
process predict_retention_times_of_possible_neoepitopes {
    publishDir "${params.outdir}/RT_prediction/"

    input:
    file txt_file_for_rt_prediction from possible_neoepitopes_list.mix(possible_neoepitopes_list_II)
    file trained_rt_param_file_II from trained_rt_params_II
    file trained_rt_set_file_II from trained_rt_set_II
    file trained_rt_model_file_II from trained_rt_model_II

    output:
    file "${txt_file_for_rt_prediction.baseName}_RTpredicted.csv" into rt_predicted_II

    when:
     params.predict_RT

    script:
    """
    RTPredict -in_text ${txt_file_for_rt_prediction} \\
              -svm_model ${trained_rt_model_file_II} \\
              -in_oligo_params ${trained_rt_param_file_II} \\
              -in_oligo_trainset ${trained_rt_set_file_II} \\
              -out_text:file "${txt_file_for_rt_prediction.baseName}_RTpredicted.csv"
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/mhcquant] Successful: $workflow.runName"
    if (!workflow.success) {
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
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

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
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/mhcquant] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nf-core/mhcquant] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/mhcquant]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/mhcquant]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/mhcquant v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
