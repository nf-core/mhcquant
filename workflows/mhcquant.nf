/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMhcquant.initialise(params, log)

// Input/outpout options
if (params.input)   { sample_sheet = file(params.input) }
if (params.fasta)   { params.fasta = params.fasta }
if (params.outdir)  { params.outdir  = './results' }

// MHC affinity prediction
if (params.predict_class_1 || params.predict_class_2)  {
    Channel.from( file(params.allele_sheet, checkIfExists: true) )
        .splitCsv( header: true, sep:'\t' )
        .multiMap { col ->
        classI: ["${col.Sample}", "${col.HLA_Alleles_Class_1}"]
        classII: ["${col.Sample}", "${col.HLA_Alleles_Class_2}"] }
        .set { ch_alleles_from_sheet }

        // Allele class 1
        if( params.predict_class_1){
            ch_alleles_from_sheet.classI
                .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
                .flatMap {it -> [tuple(it[0].toString(), it[1])] }
                .set { peptides_class_1_alleles }
        }

        // Allele class 2
        if( params.predict_class_2){
            ch_alleles_from_sheet.classII
                .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
                .flatMap { it -> [tuple(it[0].toString(), it[1])] }
                .set { peptides_class_2_alleles }
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { OPENMS_DECOYDATABASE }                                            from '../modules/local/openms_decoydatabase'
include { OPENMS_THERMORAWFILEPARSER }                                      from '../modules/local/openms_thermorawfileparser'
include { OPENMS_PEAKPICKERHIRES }                                          from '../modules/local/openms_peakpickerhires'
include { OPENMS_COMETADAPTER }                                             from '../modules/local/openms_cometadapter'
include { OPENMS_PEPTIDEINDEXER }                                           from '../modules/local/openms_peptideindexer'

include { OPENMS_FALSEDISCOVERYRATE }                                       from '../modules/local/openms_falsediscoveryrate'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_FOR_ALIGNMENT }                from '../modules/local/openms_idfilter'
include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_PSMS }                 from '../modules/local/openms_textexporter'

include { OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE }                      from '../modules/local/openms_idfilter'
include { OPENMS_IDMERGER }                                                 from '../modules/local/openms_idmerger'
include { OPENMS_PSMFEATUREEXTRACTOR }                                      from '../modules/local/openms_psmfeatureextractor'
include { OPENMS_PERCOLATORADAPTER }                                        from '../modules/local/openms_percolatoradapter'
include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_UNQUANTIFIED }         from '../modules/local/openms_textexporter'
include { CUSTOM_DUMPSOFTWAREVERSIONS }                                     from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC }                                                         from '../modules/nf-core/modules/multiqc/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK }                                                     from '../subworkflows/local/input_check'
include { INCLUDE_PROTEINS }                                                from '../subworkflows/local/include_proteins'
include { PRE_QUANTIFICATION }                                              from '../subworkflows/local/pre_quantification'
include { REFINE_FDR }                                                      from '../subworkflows/local/refine_fdr'
include { POST_QUANTIFICATION }                                             from '../subworkflows/local/post_quantification'
include { PREDICT_CLASS1 }                                                  from '../subworkflows/local/predict_class1'
include { PREDICT_CLASS2 }                                                  from '../subworkflows/local/predict_class2'
include { PREDICT_RT }                                                      from '../subworkflows/local/predict_rt'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow MHCQUANT {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Check the input file
    //
    INPUT_CHECK( params.input )
    .reads
    .set { ch_samples_from_sheet }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_samples_from_sheet
    .branch {
        meta, filename ->
            raw : meta.ext == 'raw'
                return [ meta, filename ]
            mzml :  meta.ext == 'mzml'
                return [ meta, filename ]
            other : true }
        .set { ms_files }
    // A warning message will be given when the format differs from the '.raw' or '.mzML' extention
    ms_files.other.subscribe { row -> log.warn("Unknown format for entry " + row[3] + " in provided sample sheet, line will be ignored."); exit 1 }

    // Input fasta file
    Channel.fromPath( params.fasta )
        .combine( ch_samples_from_sheet )
        .flatMap{ it -> [tuple(it[1],it[0])] }
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta }

    //
    // SUBWORKFLOW: Include protein information
    //
    if ( params.include_proteins_from_vcf ) {
        // Include the proteins from the vcf file to the fasta file
        INCLUDE_PROTEINS( input_fasta )
        ch_versions = ch_versions.mix(INCLUDE_PROTEINS.out.versions.ifEmpty(null))
        ch_fasta_file = INCLUDE_PROTEINS.out.ch_fasta_file
        ch_vcf_from_sheet = INCLUDE_PROTEINS.out.ch_vcf_from_sheet
    } else {
        ch_fasta_file = input_fasta
        ch_vcf_from_sheet = Channel.empty()
    }

    if (!params.skip_decoy_generation) {
        // Generate reversed decoy database
        OPENMS_DECOYDATABASE(ch_fasta_file)
        ch_versions = ch_versions.mix(OPENMS_DECOYDATABASE.out.versions.ifEmpty(null))
        ch_decoy_db = OPENMS_DECOYDATABASE.out.decoy
    } else {
        ch_decoy_db = ch_fasta_file
    }

    // Raw file conversion
    OPENMS_THERMORAWFILEPARSER(ms_files.raw)
    ch_versions = ch_versions.mix(OPENMS_THERMORAWFILEPARSER.out.versions.ifEmpty(null))
    if ( params.run_centroidisation ) {
        // Optional: Run Peak Picking as Preprocessing
        OPENMS_PEAKPICKERHIRES(ms_files.mzml)
        ch_versions = ch_versions.mix(OPENMS_PEAKPICKERHIRES.out.versions.ifEmpty(null))
        ch_mzml_file = OPENMS_PEAKPICKERHIRES.out.mzml
    } else {
        ch_mzml_file = ms_files.mzml
    }
    // Run comet database search
    OPENMS_COMETADAPTER(
        OPENMS_THERMORAWFILEPARSER.out.mzml
                .mix(ch_mzml_file)
                .join(ch_decoy_db, remainder:true))
    ch_versions = ch_versions.mix(OPENMS_COMETADAPTER.out.versions.ifEmpty(null))
    // Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(OPENMS_COMETADAPTER.out.idxml.join(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_PEPTIDEINDEXER.out.versions.ifEmpty(null))

    // Calculate fdr for id based alignment
    OPENMS_FALSEDISCOVERYRATE(OPENMS_PEPTIDEINDEXER.out.idxml)
    ch_versions = ch_versions.mix(OPENMS_FALSEDISCOVERYRATE.out.versions.first().ifEmpty(null))
    // Filter fdr for id based alignment
    OPENMS_IDFILTER_FOR_ALIGNMENT(OPENMS_FALSEDISCOVERYRATE.out.idxml
        .flatMap { it -> [tuple(it[0], it[1], null)]})
    ch_versions = ch_versions.mix(OPENMS_IDFILTER_FOR_ALIGNMENT.out.versions.first().ifEmpty(null))
    // Write the content to a PSMs file
    OPENMS_TEXTEXPORTER_PSMS(
        OPENMS_IDFILTER_FOR_ALIGNMENT.out.idxml
        .flatMap {
            meta, idxml ->
                ident = idxml.baseName.split('_-_')[1]
                [[[id:ident, sample:meta.sample, condition:meta.condition, ext:meta.ext], idxml]]
        }
    )

    //
    // SUBWORKFLOW: Pre-process step for the quantification of the data
    //

    if(!params.skip_quantification) {
        PRE_QUANTIFICATION(
            OPENMS_IDFILTER_FOR_ALIGNMENT.out.idxml,
            OPENMS_PEPTIDEINDEXER.out.idxml,
            ms_files.mzml
                .mix(OPENMS_THERMORAWFILEPARSER.out.mzml)
                .mix(ch_mzml_file)
        )
        ch_proceeding_idx = PRE_QUANTIFICATION.out.ch_proceeding_idx
        ch_versions = ch_versions.mix(PRE_QUANTIFICATION.out.versions.ifEmpty(null))
    } else {
        ch_proceeding_idx = OPENMS_PEPTIDEINDEXER.out.idxml
            .map {
                meta, raw ->
                [[id:meta.sample + "_" + meta.condition, sample:meta.sample, condition:meta.condition, ext:meta.ext], raw]
            }
            .groupTuple(by: [0])
    }

    // Merge aligned idXMLfiles
    OPENMS_IDMERGER(ch_proceeding_idx)
    ch_versions = ch_versions.mix(OPENMS_IDMERGER.out.versions.ifEmpty(null))
    // Extract PSM features for Percolator
    OPENMS_PSMFEATUREEXTRACTOR(OPENMS_IDMERGER.out.idxml)
    ch_versions = ch_versions.mix(OPENMS_PSMFEATUREEXTRACTOR.out.versions.ifEmpty(null))
    // Run Percolator
    OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
    ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions.ifEmpty(null))
    ch_percolator_adapter_outcome = OPENMS_PERCOLATORADAPTER.out.idxml
    // Filter by percolator q-value
    OPENMS_IDFILTER_Q_VALUE(ch_percolator_adapter_outcome.flatMap { it -> [tuple(it[0], it[1], null)]})
    ch_versions = ch_versions.mix(OPENMS_IDFILTER_Q_VALUE.out.versions.ifEmpty(null))

    //
    // SUBWORKFLOW: Refine the FDR values on the predicted subset
    //
    if ( params.refine_fdr_on_predicted_subset && params.predict_class_1 ) {
        // Run the following subworkflow
        REFINE_FDR (
            OPENMS_IDFILTER_Q_VALUE.out.idxml,
            OPENMS_PSMFEATUREEXTRACTOR.out.idxml,
            peptides_class_1_alleles
        )
        ch_versions = ch_versions.mix(REFINE_FDR.out.versions.ifEmpty(null))
        // Define the outcome of the paramer to a fixed variable
        filter_q_value = REFINE_FDR.out.filter_refined_q_value.flatMap { it -> [ tuple(it[0].sample, it[0], it[1]) ] }
    } else {
        // Make sure that the columns that consists of the ID's, sample names and the idXML file names are returned
        filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.idxml.map{ it -> [it[0].sample, it[0], it[1]] }
    }

    //
    // SUBWORKFLOW: Perform the post quantification step
    //
    if ( !params.skip_quantification) {
        POST_QUANTIFICATION (
            OPENMS_IDFILTER_FOR_ALIGNMENT.out.idxml,
            PRE_QUANTIFICATION.out.aligned_mzml,
            filter_q_value
            )
        ch_versions = ch_versions.mix(POST_QUANTIFICATION.out.versions.ifEmpty(null))
    } else {
        OPENMS_TEXTEXPORTER_UNQUANTIFIED (filter_q_value.flatMap { ident, meta, idxml -> [[meta, idxml]] })
    }

    //
    // SUBWORKFLOW: Predict class I (neoepitopes)
    //
    if ( params.predict_class_1  & !params.skip_quantification ) {
        PREDICT_CLASS1 (
            POST_QUANTIFICATION.out.mztab,
            peptides_class_1_alleles,
            ch_vcf_from_sheet
        )
        ch_versions = ch_versions.mix(PREDICT_CLASS1.out.versions.ifEmpty(null))
        ch_predicted_possible_neoepitopes = PREDICT_CLASS1.out.ch_predicted_possible_neoepitopes
    } else {
        ch_predicted_possible_neoepitopes = Channel.empty()
    }

    //
    // SUBWORKFLOW: Predict class II (neoepitopes)
    //
    if ( params.predict_class_2 & !params.skip_quantification ) {
        PREDICT_CLASS2 (
            POST_QUANTIFICATION.out.mztab,
            peptides_class_2_alleles,
            ch_vcf_from_sheet
        )
        ch_versions = ch_versions.mix(PREDICT_CLASS2.out.versions.ifEmpty(null))
        ch_predicted_possible_neoepitopes_II = PREDICT_CLASS2.out.ch_predicted_possible_neoepitopes
    } else {
        ch_predicted_possible_neoepitopes_II = Channel.empty()
    }

    //
    // SUBWORKFLOW: Predict retention time
    //
    if ( params.predict_RT ) {
        PREDICT_RT (
            filter_q_value.map{ it -> [it[1], it[2]] },
            ch_predicted_possible_neoepitopes,
            ch_predicted_possible_neoepitopes_II
        )
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowMhcquant.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        MULTIQC (
            ch_multiqc_files.collect()
        )

        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
