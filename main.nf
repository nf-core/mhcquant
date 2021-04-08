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

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    log.info nfcoreHeader()
    def command = "nextflow run nf-core/mhcquant --input 'sample_sheet.tsv' --fasta 'SWISSPROT_2020.fasta'  --allele_sheet 'alleles.tsv'  --predict_class_1  --refine_fdr_on_predicted_subset -profile standard,docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////
// Database options 
// params.include_proteins_from_vcf = false
// params.skip_decoy_generation = false


// Mass Spectronomy data processing options
// params.use_x_ions = false
// params.use_z_ions = false
// params.use_a_ions = false
// params.use_c_ions = false
// params.use_NL_ions = false
// params.remove_precursor_peak = false

// FDR scoring options
// params.refine_fdr_on_predicted_subset = false
// params.subset_affinity_threshold = 500

// Quantification options
// params.quantification_fdr = false
// params.quantification_min_prob = 0

// MHC affinity prediction options
// params.predict_class_1 = false
// params.predict_class_2 = false

// Variant options
// params.variant_reference = "GRCH38"
// params.variant_annotation_style = "SNPEFF"
// params.variant_indel_filter = false
// params.variant_frameshift_filter = false
// params.variant_snp_filter = false
////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////
// Input/outpout options
params.input = params.input ?: { log.error "Input samplesheet was not specified!"; exit 1 }()
params.outdir = params.outdir ?: { log.warn "Results into './results'.\nIf you want to define a result directory, please use the --outdir option "; return "./results" }()

// Database options
params.fasta = params.fasta ?: { log.error "No database fasta was provided, make sure you have used the '--fasta' option."; exit 1 }()
if (params.skip_decoy_generation) { 
    log.warn "Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement"
    log.warn "Decoys have to be named with DECOY_ as prefix in your fasta database" }
if (params.quantification_fdr) {
   log.warn "Quantification FDR enabled"
}

// FDR Scoring
if (params.klammer && params.description_correct_features == 0) {
    log.warn('Klammer was specified, but description of correct features was still 0. Please provide a description of correct features greater than 0.')
    log.warn('Klammer has been turned off!')
}

// Quantification options
if (params.quantification_fdr) {
   log.warn "Quantification FDR enabled"
}

////////////////////////////////////////////////////
/* --            Preprocessing                 -- */
////////////////////////////////////////////////////
sample_sheet = file(params.input)

//  Reading to content of the sample sheet and defining the different columns for further use
Channel.from( sample_sheet )
.splitCsv ( header: ['ID', 'Sample', 'Condition', 'ReplicateFileName'], sep:'\t', skip: 1) 
.map { col -> tuple("${col.ID}", "${col.Sample}", "${col.Condition}", file("${col.ReplicateFileName}", checkifExists: true))}
// .flatMap{it -> [tuple(it[0],it[1].toString(),it[2],it[3])] } 
.set { ch_samples_from_sheet }

// Checks whether the extentions of the files are known
ch_samples_from_sheet.branch {
    raw: hasExtension(it[3], 'raw')
    mzML: hasExtension(it[3], 'mzML')
    other: true
}
.set{ms_files}
// A warning message will be given when the format differs from the '.raw' or '.mzML' extention
ms_files.other.subscribe { row -> log.warn("Unknown format for entry " + row[3] + " in provided sample sheet, line will be ignored."); exit 1 }

// Define the variables depending on when the run_centroidisation parameter was called or not
if (params.run_centroidisation) {
    ms_files.mzML.set { input_mzmls_unpicked }
    Channel.empty().set { input_mzmls }
} else {
    ms_files.mzML.set { input_mzmls }
    Channel.empty().set { input_mzmls_unpicked }
}  

// Define the raw files if they are available 
ms_files.raw.set { input_raws }

// MHC affinity prediction
// TODO: test this with a allele sheet
if (params.predict_class_1 || params.predict_class_2)  {
    allele_sheet = file(params.allele_sheet, checkIfExists: true)
    
    Channel.from( allele_sheet )
        .splitCsv(header: true, sep:'\t')
        .multiMap { col -> 
        classI: ["${col.Sample}", "${col.HLA_Alleles_Class_1}"]
        classII: ["${col.Sample}", "${col.HLA_Alleles_Class_2}"] }
        .set { ch_alleles_from_sheet }
} else {
    Channel.empty().multiMap {
        sample: ""
        classI: ""
        classII: ""
    }
    .set { ch_alleles_from_sheet }
}

// Variant
if (params.include_proteins_from_vcf)  {
    vcf_sheet = file(params.vcf_sheet, checkIfExists: true)
    Channel.from( vcf_sheet )
    .splitCsv(header: ['Sample', 'VCF_FileName'], sep:'\t', skip: 1)
    .map { col -> tuple("${col.Sample}", file("${col.VCF_FileName}"),) }
    .set { ch_vcf_from_sheet }
} else {
    Channel.empty().set { ch_vcf_from_sheet }
}

////////////////////////////////////////////////////
/* --        SET CONFIGURATION VARIABLES       -- */
////////////////////////////////////////////////////

// Variant options
// variant_indel_filter = params.variant_indel_filter ? '-fINDEL' : ''
// variant_frameshift_filter = params.variant_frameshift_filter ? '-fFS' : '' 
// variant_snp_filter = params.variant_snp_filter ? '-fSNP' : ''

if (params.variant_indel_filter) { variant_indel_filter="-fINDEL" } else { variant_indel_filter="" }
if (params.variant_frameshift_filter) { variant_frameshift_filter="-fFS" } else { variant_frameshift_filter="" }
if (params.variant_snp_filter) { variant_snp_filter="-fSNP" } else { variant_snp_filter="" }


// Mass Spectronomy data processing options
x_ions = params.use_x_ions ? '-use_X_ions true' : ''
z_ions = params.use_z_ions ? '-use_Z_ions true' : ''
a_ions = params.use_a_ions ? '-use_A_ions true' : ''
c_ions = params.use_c_ions ? '-use_C_ions true' : ''
NL_ions = params.use_NL_ions ? '-use_NL_ions true' : ''
rm_precursor = params.remove_precursor_peak ? '-remove_precursor_peak true' : ''
fdr_level = (params.fdr_level == 'psm-level-fdrs') ? '' : '-'+params.fdr_level


////////////////////////////////////////////////////
/* --              CREATE CHANNELS             -- */
////////////////////////////////////////////////////
// Input fasta file

if( params.include_proteins_from_vcf) {
    Channel
        .fromPath( params.fasta )
        .combine( ch_samples_from_sheet )
        .flatMap{it -> [tuple(it[1],it[2],it[0])]}       //maps tuple to val("id"), val("Sample"), val("Condition"), val("ReplicateFileName")
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta_vcf }

    input_fasta = Channel.empty()
    input_fasta_1 = Channel.empty()
    input_fasta_2 = Channel.empty()

} else if( params.skip_decoy_generation) {
    Channel
        .fromPath( params.fasta )
        .combine( ch_samples_from_sheet )
        .flatMap{it -> [tuple(it[1],it[2],it[0])]}      //maps tuple to val("id"), val("Sample"), val("Condition"), val("ReplicateFileName")
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .into { input_fasta; input_fasta_1; input_fasta_2 }

    input_fasta_vcf = Channel.empty()

} else {
    Channel
        .fromPath( params.fasta )
        .combine( ch_samples_from_sheet )
        .flatMap{it -> [tuple(it[1].toInteger(),it[2],it[0])]}      //maps tuple to val("id"), val("Sample"), val("Condition"), val("ReplicateFileName")
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta }

    input_fasta_vcf = Channel.empty()
    input_fasta_1 = Channel.empty()
    input_fasta_2 = Channel.empty()

}
Channel
    .fromPath( params.fasta )
    .combine( ch_samples_from_sheet )
    .flatMap{it -> [tuple(it[1],it[2],it[0])] }      //maps tuple to val("id"), val("Sample"), val("Condition"), val("ReplicateFileName")
    .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
    .set { input_fasta }

// Allele class 1
if( params.predict_class_1){
    ch_alleles_from_sheet.classI
        .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
        .flatMap {it -> [tuple("id", it[0].toString(), it[1])]}     //maps tuple to val("id"), val("Sample"), val("Alleles_Class_2")
        .set { peptides_class_1_alleles }
} else {
    peptides_class_1_alleles = Channel.empty()
}

// Allele class 2
if( params.predict_class_2){

     ch_alleles_from_sheet.classII
        .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
        .flatMap {it -> [tuple("id", it[0].toString(), it[1])]}     //maps tuple to val("id"), val("Sample"), val("Alleles_Class_2")
        .set { peptides_class_2_alleles }

} else {
    peptides_class_2_alleles = Channel.empty()
}

// VCF file
if( params.include_proteins_from_vcf){
    ch_vcf_from_sheet
        .ifEmpty { exit 1, "params.vcf_sheet was empty - no vcf input file supplied" }
        .flatMap {it -> [tuple("id", it[0].toString(), it[1])]}    //maps tuple to val("id"), val("Sample"), val("VCF_FileName")
        .set { input_vcf }
} else {
    input_vcf = Channel.empty()
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

////////////////////////////////////////////////////
/* --                NF-CORE HEADER            -- */
////////////////////////////////////////////////////

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']           = custom_runName ?: workflow.runName
summary['Pipeline Name']      = 'nf-core/mhcquant'
summary['Pipeline Version']   = workflow.manifest.version
summary['Run Name']           = custom_runName ?: workflow.runName
summary['MS Samples']         = input_mzmls.mix(input_raws)
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
summary['Launch dir'] = workflow.launchDir
summary['User'] = workflow.userName
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
//////////////////////////////////////////////////////////////////////////////////////////////
/*
 * STEP 0 - Output Description HTML
*/ 
// process output_documentation {
//     publishDir "${params.outdir}/pipeline_info", mode: 'copy'

//     input:
//     file ch_output_docs //from ch_output_docs

//     output:
//     file "results_description.html"

//     script:
//     """
//     markdown_to_html.py $ch_output_docs -o results_description.html
//     ""

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

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def modules = params.modules.clone()
def openms_map_aligner_identification_options = modules['openms_map_aligner_identification']
def openms_comet_adapter_options = modules['openms_comet_adapter']
def generate_proteins_from_vcf_options = modules['generate_proteins_from_vcf']
// def quantify_id_targeted_options = modules['quantify_id_targeted']

openms_comet_adapter_options.args += x_ions + z_ions + c_ions + a_ions + NL_ions + rm_precursor
generate_proteins_from_vcf_options.args += variant_indel_filter + variant_snp_filter + variant_frameshift_filter

// include { GET_SOFTWARE_VERSIONS }                           from './modules/local/process/get_software_versions'                            addParams( options: [publish_files : ['tsv':'']]   )
include { INPUT_CHECK }                                     from './modules/local/subworkflow/input_check'                                  addParams( options: [:] )           
include { GENERATE_PROTEINS_FROM_VCF }                      from './modules/local/process/generate_proteins_from_vcf'                       addParams( options: generate_proteins_from_vcf_options )
include { OPENMS_DECOYDATABASE }                            from './modules/local/process/openms_decoydatabase'                             addParams( options: [:] )
include { OPENMS_THERMORAWFILEPARSER }                      from './modules/local/process/openms_thermorawfileparser'                       addParams( options: [:] )
include { OPENMS_PEAKPICKERHIRES }                          from './modules/local/process/openms_peakpickerhires'                           addParams( options: [:] )
include { OPENMS_COMETADAPTER }                             from './modules/local/process/openms_cometadapter'                              addParams( options: openms_comet_adapter_options )
include { OPENMS_PEPTIDEINDEXER }                           from './modules/local/process/openms_peptideindexer'                            addParams( options: [:] )
include { OPENMS_FALSEDISCOVERYRATE }                       from './modules/local/process/openms_falsediscoveryrate'                        addParams( options: [:] )
include { FILTER_FDR_FOR_ID_ALIGNMENT }                     from './modules/local/process/filter_fdr_for_id_alignment'                      addParams( options: [:] )
include { OPENMS_MAPALIGNERIDENTIFICATION }                 from './modules/local/process/openms_mapaligneridentification'                  addParams( options: openms_map_aligner_identification_options )

include { 
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERMZML
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERIDXML }        from './modules/local/process/openms_maprttransformer'                 addParams( options: [:] )

include { OPENMS_IDMERGER }                                 from './modules/local/process/openms_idmerger'                                  addParams( options: [:] )
include { OPENMS_PSMFEATUREEXTRACTOR }                      from './modules/local/process/openms_psmfeatureextractor'                       addParams( options: [:] )
include { RUN_PERCOLATOR }                                  from './modules/local/process/run_percolator'                                   addParams( options: [:] )
include { FILTER_BY_Q_VALUE }                               from './modules/local/process/filter_by_q_value'                                addParams( options: [:] )

include { REFINE_FDR_ON_PREDICTED_SUBSET }                  from './modules/local/subworkflow/refine_fdr_on_predicted_subset'               addParams( options: [:] )

include { OPENMS_FEATUREFINDERIDENTIFICATION }              from './modules/local/process/openms_featurefinderidentification'               addParams( options: [:] )
include { OPENMS_FEATURELINKERUNLABELEDKD }                 from './modules/local/process/openms_featurelinkerunlabeledkd'                  addParams( options: [:] )
include { OPENMS_IDCONFLICTRESOLVER }                       from './modules/local/process/openms_idconflictresolver'                        addParams( options: [:] )
include { OPENMS_TEXTEXPORTER }                             from './modules/local/process/openms_textexporter'                              addParams( options: [:] )
include { OPENMS_MZTABEXPORTER }                            from './modules/local/process/openms_mztabexporter'                             addParams( options: [:] )

include { PREDICT_PEPTIDES_MHCFLURRY_CLASS_1 }              from './modules/local/process/predict_peptides_mhcflurry_class_1'               addParams( options: [:] )
include { PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 }          from './modules/local/process/preprocess_peptides_mhcnuggets_class_2'           addParams( options: [:] )
include { PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2 }             from './modules/local/process/predict_peptides_mhcnuggets_class_2'              addParams( options: [:] )
include { POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 }         from './modules/local/process/postprocess_peptides_mhcnuggets_class_2'          addParams( options: [:] )
include { PREDICT_POSSIBLE_NEOEPITOPES }                    from './modules/local/process/predict_possible_neoepitopes'                     addParams( options: [:] )
include { PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES }            from './modules/local/process/predict_possible_class_2_neoepitopes'             addParams( options: [:] )
include { RESOLVE_FOUND_NEOEPITOPES }                       from './modules/local/process/resolve_found_neoepitopes'                        addParams( options: [:] )
include { RESOLVE_FOUND_CLASS_2_NEOEPITOPES }               from './modules/local/process/resolve_found_class_2_neoepitopes'                addParams( options: [:] )
include { PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1 }           from './modules/local/process/predict_neoepitopes_mhcflurry_class_1'            addParams( options: [:] )
include { PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 }       from './modules/local/process/preprocess_neoepitopes_mhcnuggets_class_2'        addParams( options: [:] )
include { PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2 }          from './modules/local/process/predict_neoepitopes_mhcnuggets_class_2'           addParams( options: [:] )
include { POSTPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 }      from './modules/local/process/postprocess_neoepitopes_mhcnuggets_class_2'       addParams( options: [:] )
include { OPENMS_RTMODEL }                                  from './modules/local/process/openms_rtmodel'                                   addParams( options: [:] )
include { PREDICT_RETENTION_TIMES_OF_FOUND_PEPTIDES }       from './modules/local/process/predict_retention_times_of_found_peptides'        addParams( options: [:] )
include { PREDICT_RETENTION_TIMES_OF_POSSIBLE_NEOEPITOPES } from './modules/local/process/predict_retention_times_of_possible_neoepitopes'  addParams( options: [:] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow {
    // TODO: start with the generation of the subworkflows to call the beloning tool versions 
    //////////////////////////////////////////////////////////////////////////////////////////////
    // TODO: reduce the number of modules (combine the same functions into one)
    //////////////////////////////////////////////////////////////////////////////////////////////
    ch_software_versions = Channel.empty()
    // // Example: ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))
    // GET_SOFTWARE_VERSIONS (ch_software_versions.map { it }.collect())

    // TODO: Include properly
    // DONE: Input based on the meta id
    // INPUT_CHECK(sample_sheet)
        // .map { meta, raw -> meta.id = meta.id.split('_')[0..-1].join('_')
        // [ meta, raw ] }
        // .groupTuple(by: [0])
        // .map { it -> [ it[0], it[1].flatten() ] }
        // .set { ch_samples_from_sheet }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // If specified translate variants to proteins and include in reference fasta
    GENERATE_PROTEINS_FROM_VCF(input_fasta.combine(input_vcf, by:1))
    ch_software_versions = ch_software_versions.mix(GENERATE_PROTEINS_FROM_VCF.out[1].ifEmpty(null))
    //Generate reversed decoy database
    OPENMS_DECOYDATABASE(input_fasta.mix(GENERATE_PROTEINS_FROM_VCF.out[0]))
    ch_software_versions = ch_software_versions.mix(OPENMS_DECOYDATABASE.out[1].first().ifEmpty(null))
    // Raw file conversion
    OPENMS_THERMORAWFILEPARSER(input_raws)
    ch_software_versions = ch_software_versions.mix(OPENMS_THERMORAWFILEPARSER.out[1].first().ifEmpty(null))
    ////////////////////////////////////////////////////////////////////////////////////////////
    //TODO: create if statement, remove when from file
    //Optional: Run Peak Picking as Preprocessing
    OPENMS_PEAKPICKERHIRES(input_mzmls_unpicked)
    ch_software_versions = ch_software_versions.mix(OPENMS_PEAKPICKERHIRES.out[1].first().ifEmpty(null))
    //Run comet database search
    OPENMS_COMETADAPTER(OPENMS_THERMORAWFILEPARSER.out[0].mix(input_mzmls.mix(OPENMS_PEAKPICKERHIRES.out[0])).join(OPENMS_DECOYDATABASE.out[0].mix(input_fasta_1), by:1, remainder:true))
    ch_software_versions = ch_software_versions.mix(OPENMS_COMETADAPTER.out[1].first().ifEmpty(null))
    //Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(OPENMS_COMETADAPTER.out[0].join(OPENMS_DECOYDATABASE.out[0].mix(input_fasta_2), by:1))
    ch_software_versions = ch_software_versions.mix(OPENMS_PEPTIDEINDEXER.out[1].first().ifEmpty(null))
    //////////////////////////////////////////////////////////////////////////////////////////
     //Calculate fdr for id based alignment
     //TODO: create if statement, remove when from file
    OPENMS_FALSEDISCOVERYRATE(OPENMS_PEPTIDEINDEXER.out[0])
    ch_software_versions = ch_software_versions.mix(OPENMS_FALSEDISCOVERYRATE.out[1].first().ifEmpty(null))
     //Filter fdr for id based alignment
    FILTER_FDR_FOR_ID_ALIGNMENT(OPENMS_FALSEDISCOVERYRATE.out.idxml) 
    ch_software_versions = ch_software_versions.mix(FILTER_FDR_FOR_ID_ALIGNMENT.out[1].first().ifEmpty(null))
     //Compute alignment rt transformation
    OPENMS_MAPALIGNERIDENTIFICATION(FILTER_FDR_FOR_ID_ALIGNMENT.out.idxml)
    ch_software_versions = ch_software_versions.mix(OPENMS_MAPALIGNERIDENTIFICATION.out[1].first().ifEmpty(null))

    //Intermediate Step to join RT transformation files with mzml and idxml channels
    if(!params.skip_quantification) {
        input_mzmls
        .mix(OPENMS_THERMORAWFILEPARSER.out.mzml)
        .mix(OPENMS_PEAKPICKERHIRES.out.mzml)
        .flatMap { it -> [tuple(it[0].toInteger(), it[1], it[2], it[3])]}
        .join(OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml.transpose().flatMap{ it -> [tuple(it[1].baseName.split('_-_')[0].toInteger(), it[0], it[1])]}, by: [0,1])
        .set{joined_trafos_mzmls}

        OPENMS_PEPTIDEINDEXER.out.idxml
        .flatMap { it -> [tuple(it[0].toInteger(), it[1], it[2], it[3])]}
        .join(OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml.transpose().flatMap{ it -> [tuple(it[1].baseName.split('_-_')[0].toInteger(), it[0], it[1])]}, by: [0,1])
        .set{joined_trafos_ids}

        id_files_idx_original = Channel.empty()
    } else {
        joined_trafos_mzmls = Channel.empty()
        joined_trafos_ids = Channel.empty()
        id_files_idx_original = OPENMS_PEPTIDEINDEXER.out.idxml
    }

    //Align mzML files using trafoXMLs
    OPENMS_MAPRTTRANSFORMERMZML(joined_trafos_mzmls)
    ch_software_versions = ch_software_versions.mix(OPENMS_MAPRTTRANSFORMERMZML.out[1].first().ifEmpty(null))
    //Align unfiltered idXMLfiles using trafoXMLs
    OPENMS_MAPRTTRANSFORMERIDXML(joined_trafos_ids)
    ch_software_versions = ch_software_versions.mix(OPENMS_MAPRTTRANSFORMERIDXML.out[1].first().ifEmpty(null))
    ////////////////////////////////////////////////////////////////////////////////////////////
    //Merge aligned idXMLfiles
    OPENMS_IDMERGER(OPENMS_MAPRTTRANSFORMERIDXML.out[0].mix(id_files_idx_original).groupTuple(by: 1))
    ch_software_versions = ch_software_versions.mix(OPENMS_IDMERGER.out[1].first().ifEmpty(null))
    //Extract PSM features for Percolator
    OPENMS_PSMFEATUREEXTRACTOR(OPENMS_IDMERGER.out.idxml)
    ch_software_versions = ch_software_versions.mix(OPENMS_PSMFEATUREEXTRACTOR.out[1].first().ifEmpty(null))
    //Run Percolator
    RUN_PERCOLATOR(OPENMS_PSMFEATUREEXTRACTOR.out.idxml, fdr_level)
    ch_software_versions = ch_software_versions.mix(RUN_PERCOLATOR.out[1].first().ifEmpty(null))
    //Filter by percolator q-value
    FILTER_BY_Q_VALUE(RUN_PERCOLATOR.out.idxml) // NOTE: Same function was used for (!)params.refine_fdr_on_predicted_subset
    ch_software_versions = ch_software_versions.mix(FILTER_BY_Q_VALUE.out[1].first().ifEmpty(null))
    ///////////////////////////////////////////////////////////////////////////////////////
    //Refine_fdr_on_predicted_subset
    //TODO: change the different parameters to an options.args
    if ( params.refine_fdr_on_predicted_subset) {
       // Run the following subworkflow
        REFINE_FDR_ON_PREDICTED_SUBSET (
            FILTER_BY_Q_VALUE.out[0], // Added the condition column... might become obsolete later
            OPENMS_PSMFEATUREEXTRACTOR.out[0],
            peptides_class_1_alleles,
            fdr_level // TODO: Add the following parameters to a option.args
        )
       // Define the outcome of the paramer to a fixed variable
        // filter_refined_q_value = REFINE_FDR_ON_PREDICTED_SUBSET.out[0]  
    // } else {
        // filter_refined_q_value = Channel.empty()
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////  
    //FILTER_FDR_FOR_ID_ALIGNMENT.out[0]
    //.flatMap { it -> [ tuple( it[0], it[1], it[2], it[3] ) ] }
    //.join( OPENMS_MAPRTTRANSFORMERMZML.out[0], by: [0,1,2] )
    //.combine( FILTER_BY_Q_VALUE.out[0].flatMap { it -> [ tuple( it[0], it[1], it[3] ) ] }
    //    .mix( filter_refined_q_value ), by: 1 )
    //.set{ joined_mzmls_ids_quant }
    //// .ifEmpty('Hello')
//
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //if (! params.skip_quantification) {
    //    // Quantify identifications using targeted feature extraction
    //    OPENMS_FEATUREFINDERIDENTIFICATION(joined_mzmls_ids_quant)
    //    // Link extracted features
    //    OPENMS_FEATURELINKERUNLABELEDKD(OPENMS_FEATUREFINDERIDENTIFICATION.out.featurexml.groupTuple(by:1))  
    //    // Resolve conflicting ids matching to the same feature
    //    OPENMS_IDCONFLICTRESOLVER(OPENMS_FEATURELINKERUNLABELEDKD.out.consensusxml)
    //}
    //////////////////////////////////////////////////////////////////////////////////////////////////
//
    //// Export all information as text to csv
    //OPENMS_TEXTEXPORTER(OPENMS_IDCONFLICTRESOLVER.out.consensusxml)
    //// Export all information as mzTab
    //OPENMS_MZTABEXPORTER(
    //    OPENMS_IDCONFLICTRESOLVER.out.consensusxml
    //    .flatMap { it -> [tuple("id", it[0], it[0], it[1])]}
    //)



     //////////////////////////////////////////////////////////////////////////////////////////////

    ////  TODO: Make a subworkflow, proceed with the new container
    ////  START:
    //// If specified predict peptides using MHCFlurry
    //PREDICT_PEPTIDES_MHCFLURRY_CLASS_1(OPENMS_MZTABEXPORTER.out.mztab.combine(peptides_class_1_alleles, by:1)) // Note for replacement epytope
    //// Preprocess found peptides for MHCNuggets prediction class 2
    //PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2(OPENMS_MZTABEXPORTER.out.mztab) // Note for replacement epytope
    //// Predict found peptides using MHCNuggets class 2
    //PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2(PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out[0].join(peptides_class_2_alleles, by:1)) // Note for replacement epytope
    //// Postprocess predicted MHCNuggets peptides class 2
    //POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2(PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2.out.csv.join(PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out[1], by:1))  // Note for replacement epytope
//
    //// Predict all possible neoepitopes from vcf
    //PREDICT_POSSIBLE_NEOEPITOPES(peptides_class_1_alleles.join(input_vcf, by:[0,1], remainder:true)) // Note for replacement epytope
    //// Predict all possible class 2 neoepitopes from vcf
    //PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES(peptides_class_2_alleles.join(input_vcf, by:[0,1], remainder:true)) // Note for replacement epytope
    //// Resolve found neoepitopes
    //RESOLVE_FOUND_NEOEPITOPES(OPENMS_MZTABEXPORTER.out.mztab.join(PREDICT_POSSIBLE_NEOEPITOPES.out[0], by:[0,1], remainder:true))
    //// Resolve found class 2 neoepitopes
    //RESOLVE_FOUND_CLASS_2_NEOEPITOPES(OPENMS_MZTABEXPORTER.out.mztab.join(PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES.out[0], by:[0,1], remainder:true)) 
    //// Predict class 1 neoepitopes MHCFlurry
    //PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1(peptides_class_1_alleles.join(RESOLVE_FOUND_NEOEPITOPES.out.csv, by:1)) // Note for replacement epytope
    //// Preprocess resolved neoepitopes in a format that MHCNuggets understands
    //PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2(RESOLVE_FOUND_CLASS_2_NEOEPITOPES.out.csv) // Note for replacement epytope
    //// Predict class 2 MHCNuggets
    //PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2(PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.preprocessed.join(peptides_class_2_alleles, by:1)) // Note for replacement epytope
    //// Class 2 MHCNuggets Postprocessing
    //POSTPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2(RESOLVE_FOUND_CLASS_2_NEOEPITOPES.out.csv.join(PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.csv, by:1)) // Note for replacement epytope
//
    //// Train Retention Times Predictor
    //OPENMS_RTMODEL(FILTER_BY_Q_VALUE.out.idxml.mix(REFINE_FDR_ON_PREDICTED_SUBSET.out.filter_refined_q_value))
    //// Retention Times Predictor Found Peptides
    //PREDICT_RETENTION_TIMES_OF_FOUND_PEPTIDES(FILTER_BY_Q_VALUE.out.idxml.mix(REFINE_FDR_ON_PREDICTED_SUBSET.out.filter_refined_q_value).join(OPENMS_RTMODEL.out.complete, by:[0,1]))  
    //// Retention Times Predictor possible Neoepitopes
    //PREDICT_RETENTION_TIMES_OF_POSSIBLE_NEOEPITOPES(PREDICT_POSSIBLE_NEOEPITOPES.out[1].mix(PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES.out[1]).join(OPENMS_RTMODEL.out.complete, by:[0,1])) 
    //// END
//

}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////
// workflow.onComplete {
    // Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
    // Completion.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
// }
// 
////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////