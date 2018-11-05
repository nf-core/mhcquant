#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/openmspeptidequant
========================================================================================
 nf-core/openmspeptidequant Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/openmspeptidequant
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/openmspeptidequant v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/openmspeptidequant --mzmls '*.mzML' --fasta '*.fasta' -profile standard,docker

    Mandatory arguments:
      --mzmls                       Path to input data (must be surrounded with quotes)
      --fasta                       Path to Fasta reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:
      --precursor_mass_tolerance    Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance     Mass tolerance of fragment mass bin (ppm)
      --fragment_bin_offset         Offset of fragment mass bin (Comet specific parameter)
      --fdr_threshold               Threshold for FDR filtering
      --fdr_level                   Level of FDR calculation ('peptide-level-fdrs', 'psm-level-fdrs', 'protein-level-fdrs')
      --digest_mass_range           Mass range of peptides considered for matching
      --activation_method           Fragmentation method ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')
      --enzyme                      Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see OpenMS enzymes)
      --number_mods                 Maximum number of modifications of PSMs
      --fixed_mods                  Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods               Variable modifications ('Oxidation (M)', see OpenMS modifications)

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
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
params.mzmls ?: params.mzmlPaths ?: { log.error "No read data privided. Make sure you have used the '--mzmls' option."; exit 1 }()
params.fasta ?: params.fastaPath ?: { log.error "No read data privided. Make sure you have used the '--fasta' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()




/ Parameters
params.num_threads = 5

params.fmt = 0.02
params.pmt = 5
params.fbo = 0
params.fdr = 0.01
params.maxmod = 3

params.num_hits = 1
params.dmr = "800:2500"
params.msLevels = 2
params.centroided = "False"

params.prec_charge = '2:3'
params.activ_method = 'ALL'

params.enzyme = 'unspecific cleavage'
params.fixed_mods = ''
params.variable_mods = 'Oxidation (M)'


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

//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


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
 * Create a channel for input read files
 */
 if(params.readPaths){
     if(params.singleEnd){
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_trimming }
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_trimming }
     }
 } else {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming }
 }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/openmspeptidequant v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/openmspeptidequant'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
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
    id: 'nf-core-openmspeptidequant-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/openmspeptidequant Workflow Summary'
    section_href: 'https://github.com/nf-core/openmspeptidequant'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - generate reversed decoy database
 */
process generate_decoy_database {
    publishDir "${results_path}/"

    input:
     file fastafile
    output:
     file "${fastafile.baseName}_decoy.fasta" into fastafile_decoy
     
    script:
     """
     ${openms_path}/DecoyDatabase  -in ${fastafile} -out ${fastafile.baseName}_decoy.fasta -decoy_string DECOY_ -decoy_string_position prefix
     """
}


/*
 * STEP 2 - run comet database search
 */
process db_search_comet {
    publishDir "${results_path}/"
 
    input:
     set mzmlID, file(mzml_file) from mzmlfiles_search
     file fasta_decoy from fastafile_decoy

    output:
     set mzmlID, file("${mzmlID}.idXML") into id_files

    script:
     """
     ${openms_path}/CometAdapter -comet_executable ${comet_executable_path} -in ${mzml_file} -out ${mzmlID}.idXML -threads ${params.num_threads} -database ${fasta_decoy} -precursor_mass_tolerance ${params.pmt} -fragment_bin_tolerance ${params.fmt} -fragment_bin_offset ${params.fbo} -num_hits ${params.num_hits} -digest_mass_range ${params.dmr} -max_variable_mods_in_peptide ${params.maxmod} -allowed_missed_cleavages 0 -precursor_charge ${params.prec_charge} -activation_method ${params.activ_method} -use_NL_ions true -variable_modifications '${params.variable_mods}' -fixed_modifications ${params.fixed_mods} -enzyme '${params.enzyme}'
     """

}


/*
 * STEP 3 - index decoy and target hits
 */
process index_peptides {
    publishDir "${results_path}/"
 
    input:
     set ID_idx, file(id_file) from id_files
     file fasta_decoy from fastafile_decoy

    output:
     set ID_idx, file("${ID_idx}_idx.idXML") into id_files_idx, id_files_idx_original

    script:
     """
     ${openms_path}/PeptideIndexer -in ${id_file} -out ${ID_idx}_idx.idXML -threads ${params.num_threads} -fasta ${fasta_decoy} -decoy_string DECOY -enzyme:specificity none
     """

}


/*
 * STEP 4 - calculate fdr for id based alignment
 */
process calculate_fdr_for_idalignment {
    publishDir "${results_path}/"
 
    input:
     set ID_idx_fdr, file(id_file_idx) from id_files_idx

    output:
     set ID_idx_fdr, file("${ID_idx_fdr}_fdr.idXML") into id_files_idx_fdr

    script:
     """
     ${openms_path}/FalseDiscoveryRate -in ${id_file_idx} -out ${ID_idx_fdr}_fdr.idXML -threads ${params.num_threads}
     """

}


/*
 * STEP 5 - filter fdr for id based alignment
 */
process filter_fdr_for_idalignment {
    publishDir "${results_path}/"
 
    input:
     set ID_idx_fdr_filtered, file(id_file_idx_fdr) from id_files_idx_fdr

    output:
     set ID_idx_fdr_filtered, file("${ID_idx_fdr_filtered}_filtered.idXML") into id_files_idx_fdr_filtered

    script:
     """
     ${openms_path}/IDFilter -in ${id_file_idx_fdr} -out ${ID_idx_fdr_filtered}_filtered.idXML -threads ${params.num_threads} -score:pep 0.05  -remove_decoys
     """

}


/*
 * STEP 6 - compute alignment rt transformation
 */
process align_ids {
   publishDir "${results_path}/"

    input:
     file id_names from id_files_idx_fdr_filtered.collect{it[1]}

    output:
     file '*.trafoXML' into id_files_trafo_mzml, id_files_trafo_idxml

    script:
     def out_names = id_names.collect { it.baseName+'.trafoXML' }.join(' ')
     """
     ${openms_path}/MapAlignerIdentification -in $id_names -trafo_out $out_names
     """

}


/*
 * STEP 7 - align mzML files using trafoXMLs
 */
process align_mzml_files {
    publishDir "${results_path}/"

    input:
     set ID_trafo_mzml, file(id_file_trafo) from id_files_trafo_mzml
     set mzmlID_align, file(mzml_file_align) from mzmlfiles_align

    output:
     set mzmlID_align, file("${mzmlID_align}_aligned.mzML") into mzml_files_aligned

    script:
     """
     ${openms_path}/MapRTTransformer -in ${mzml_file_align} -trafo_in ${id_file_trafo} -out ${mzmlID_align}_aligned.mzML -threads ${params.num_threads}
     """

}


/*
 * STEP 8 - align unfiltered idXMLfiles using trafoXMLs
 */
process align_idxml_files {
    publishDir "${results_path}/"

    input:
     set ID_trafo_idxml, file(idxml_file_trafo) from id_files_trafo_idxml
     set ID_align, file(idxml_file_align) from id_files_idx_original

    output:
     set ID_align, file("${ID_align}_aligned.idXML") into idxml_files_aligned

    script:
     """
     ${openms_path}/MapRTTransformer -in ${idxml_file_align} -trafo_in ${idxml_file_trafo} -out ${ID_align}_aligned.idXML -threads ${params.num_threads}
     """

}


/*
 * STEP 9 - merge aligned idXMLfiles
 */
process merge_aligned_idxml_files {
    publishDir "${results_path}/"

    input:
     //stdin idxml_files_aligned.collect()
     file ids_aligned from idxml_files_aligned.flatMap().buffer( size: 2 ).collect{ it[1] }

    output:
     file "all_ids_merged.idXML" into id_merged
    
    script:
     """
     ${openms_path}/IDMerger -in $ids_aligned -out all_ids_merged.idXML -threads ${params.num_threads}  -annotate_file_origin
     """

}


/*
 * STEP 10 - extract PSM features for Percolator
 */
process extract_psm_features_for_percolator {
    publishDir "${results_path}/"
 
    input:
     set ID_psm, file(id_file_merged) from id_merged.map { file -> tuple(file.baseName, file)}

    output:
     set ID_psm, file("${ID_psm}_psm.idXML") into id_files_merged_psm

    script:
     """
     ${openms_path}/PSMFeatureExtractor -in ${id_file_merged} -out ${ID_psm}_psm.idXML -threads ${params.num_threads} 
     """

}


/*
 * STEP 11 - run Percolator
 */
process run_percolator {
    publishDir "${results_path}/"
 
    input:
     set ID_perc, file(id_file_psm) from id_files_merged_psm

    output:
     set ID_perc, file("${ID_perc}_psm_perc.idXML") into id_files_merged_psm_perc

    script:
     """
     ${openms_path}/PercolatorAdapter -in ${id_file_psm} -out ${ID_perc}_psm_perc.idXML -threads ${params.num_threads} -enzyme no_enzyme -percolator_executable ${percolator_executable_path}
     """

}


/*
 * STEP 12 - filter by percolator q-value
 */
process filter_q_value {
    publishDir "${results_path}/"
 
    input:
     set ID_perc_filtered, file(id_file_perc) from id_files_merged_psm_perc

    output:
     set ID_perc_filtered, file("${ID_perc_filtered}_psm_perc_filtered.idXML") into id_files_merged_psm_perc_filtered

    script:
     """
     ${openms_path}/IDFilter -in ${id_file_perc} -out ${ID_perc_filtered}_psm_perc_filtered.idXML -threads ${params.num_threads} -score:pep 9999  -remove_decoys
     """

}


/*
 * STEP 13 - quantify identifications using targeted feature extraction
 */
process quantify_identifications_targeted {
    publishDir "${results_path}/"
 
    input:
     set ID_quant, file(id_file_quant) from id_files_merged_psm_perc_filtered.first()
     set FEAT_quant, file(mzml_quant) from mzml_files_aligned

    output:
     set FEAT_quant, file("${FEAT_quant}.featureXML") into feature_files

    script:
     """
     ${openms_path}/FeatureFinderIdentification -in ${mzml_quant} -id ${id_file_quant} -out ${FEAT_quant}.featureXML -threads ${params.num_threads}
     """

}


/*
 * STEP 14 - link extracted features
 */
process link_extracted_features {
    publishDir "${results_path}/"

    input:
     file feautres from feature_files.collect{it[1]}

    output:
     file "all_features_merged.consensusXML" into consensus_file
    
    script:
     """
     ${openms_path}/FeatureLinkerUnlabeledKD -in $feautres -out 'all_features_merged.consensusXML' -threads ${params.num_threads}
     """

}


/*
 * STEP 15 - resolve conflicting ids matching to the same feature
 */
process resolve_conflicts {
    publishDir "${results_path}/"
 
    input:
     set CONS, file(consensus) from consensus_file.map { file -> tuple(file.baseName, file)}

    output:
     set CONS, file("${CONS}_resolved.consensusXML") into consensus_file_resolved

    script:
     """
     ${openms_path}/IDConflictResolver -in ${consensus} -out ${CONS}_resolved.consensusXML -threads ${params.num_threads}
     """

}


/*
 * STEP 16 - export all information as text to csv
 */
process export_text {
    publishDir "${results_path}/"
 
    input:
     set CONS_resolved, file(consensus_resolved) from consensus_file_resolved

    output:
     set CONS_resolved, file("${CONS_resolved}.csv") into consensus_text

    script:
     """
     ${openms_path}/TextExporter -in ${consensus_resolved} -out ${CONS_resolved}.csv -threads ${params.num_threads} -id:add_hit_metavalues 0 -id:add_metavalues 0 -id:peptides_only
     """

}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/openmspeptidequant] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/openmspeptidequant] FAILED: $workflow.runName"
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
          log.info "[nf-core/openmspeptidequant] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/openmspeptidequant] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/openmspeptidequant] Pipeline Complete"

}
