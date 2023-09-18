/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMhcquant.initialise(params, log)

// Input/output options
if (params.input)   { sample_sheet = file(params.input) }
if (params.fasta)   { params.fasta = params.fasta }

// MHC affinity prediction
if (params.predict_class_1 || params.predict_class_2) {
    Channel.from(file(params.allele_sheet, checkIfExists: true))
        .splitCsv(header: true, sep:'\t')
        .multiMap { col ->
            classI: ["${col.Sample}", "${col.HLA_Alleles_Class_1}"]
            classII: ["${col.Sample}", "${col.HLA_Alleles_Class_2}"] }
        .set { ch_alleles_from_sheet }

        // Allele class 1
        if (params.predict_class_1) {
            ch_alleles_from_sheet.classI
                .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
                .flatMap { it -> [tuple(it[0].toString(), it[1])] }
                .set { peptides_class_1_alleles }
        }

        // Allele class 2
        if (params.predict_class_2) {
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

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { OPENMS_DECOYDATABASE }                                            from '../modules/local/openms_decoydatabase'
include { THERMORAWFILEPARSER }                                             from '../modules/local/thermorawfileparser'
include { TDF2MZML }                                                        from '../modules/local/tdf2mzml'
include { OPENMS_PEAKPICKERHIRES }                                          from '../modules/local/openms_peakpickerhires'
include { OPENMS_FILEFILTER }                                               from '../modules/local/openms_filefilter'
include { OPENMS_COMETADAPTER }                                             from '../modules/local/openms_cometadapter'
include { OPENMS_PEPTIDEINDEXER }                                           from '../modules/local/openms_peptideindexer'
include { DEEPLC }                                                          from '../modules/local/deeplc'
include { MS2PIP }                                                          from '../modules/local/ms2pip'

include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_COMET }                from '../modules/local/openms_textexporter'

include { OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE }                      from '../modules/local/openms_idfilter'
include { OPENMS_IDMERGER }                                                 from '../modules/local/openms_idmerger'
include { OPENMS_IDMERGER as OPENMS_IDMERGER_QUANT }                        from '../modules/local/openms_idmerger'
include { OPENMS_IDRIPPER }                                                 from '../modules/local/openms_idripper'

include { OPENMS_PSMFEATUREEXTRACTOR }                                      from '../modules/local/openms_psmfeatureextractor'
include { OPENMS_PERCOLATORADAPTER }                                        from '../modules/local/openms_percolatoradapter'
include { PYOPENMS_IDFILTER }                                               from '../modules/local/pyopenms_idfilter'
include { PYOPENMS_IONANNOTATOR }                                           from '../modules/local/pyopenms_ionannotator'

include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_ID }                  from '../modules/local/openms_textexporter'
include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_QUANTIFIED }          from '../modules/local/openms_textexporter'
include { OPENMS_MZTABEXPORTER as OPENMS_MZTABEXPORTER_ID }                from '../modules/local/openms_mztabexporter'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

include { INCLUDE_PROTEINS }                                                from '../subworkflows/local/include_proteins'
include { MAP_ALIGNMENT }                                                   from '../subworkflows/local/map_alignment'
include { REFINE_FDR }                                                      from '../subworkflows/local/refine_fdr'
include { PROCESS_FEATURE }                                                 from '../subworkflows/local/process_feature.nf'
include { PREDICT_CLASS1 }                                                  from '../subworkflows/local/predict_class1'
include { PREDICT_CLASS2 }                                                  from '../subworkflows/local/predict_class2'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow MHCQUANT {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Check the input file
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    INPUT_CHECK.out.ms_runs
        .branch {
            meta, filename ->
                raw : meta.ext == 'raw'
                    return [ meta, filename ]
                mzml : meta.ext == 'mzml'
                    return [ meta, filename ]
                tdf : meta.ext == 'd'
                    return [ meta, filename ]
                other : true }
        .set { branched_ms_files }

    // Input fasta file
    Channel.fromPath(params.fasta)
        .combine(INPUT_CHECK.out.ms_runs)
        .flatMap{ it -> [tuple(it[1],it[0])] }
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta }

    //
    // SUBWORKFLOW: Include protein information
    //
    if (params.include_proteins_from_vcf) {
        // Include the proteins from the vcf file to the fasta file
        INCLUDE_PROTEINS(input_fasta)
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

    ch_ms_files = (branched_ms_files.mzml)
    // Raw file conversion
    THERMORAWFILEPARSER(branched_ms_files.raw)
    ch_versions = ch_versions.mix(THERMORAWFILEPARSER.out.versions.ifEmpty(null))
    ch_ms_files = ch_ms_files.mix(THERMORAWFILEPARSER.out.mzml)

    // timsTOF data conversion
    TDF2MZML(branched_ms_files.tdf)
    ch_versions = ch_versions.mix(TDF2MZML.out.versions.ifEmpty(null))
    ch_ms_files = ch_ms_files.mix(TDF2MZML.out.mzml)

    if (params.run_centroidisation) {
        // Optional: Run Peak Picking as Preprocessing
        OPENMS_PEAKPICKERHIRES(ch_ms_files)
        ch_versions = ch_versions.mix(OPENMS_PEAKPICKERHIRES.out.versions.ifEmpty(null))
        ch_mzml_file = OPENMS_PEAKPICKERHIRES.out.mzml
    } else {
        ch_mzml_file = ch_ms_files
    }

    // Clean up mzML files
    OPENMS_FILEFILTER(ch_mzml_file)
    ch_versions = ch_versions.mix(OPENMS_FILEFILTER.out.versions.ifEmpty(null))
    ch_clean_mzml_file = OPENMS_FILEFILTER.out.cleaned_mzml

    // Run comet database search
    OPENMS_COMETADAPTER(
            ch_clean_mzml_file.join(ch_decoy_db, remainder:true))

    // Run DeepLC if specified
    if (params.use_deeplc){
        DEEPLC(OPENMS_COMETADAPTER.out.idxml)
        ch_versions = ch_versions.mix(DEEPLC.out.versions.ifEmpty(null))
        ch_comet_out_idxml = DEEPLC.out.idxml
    } else {
        ch_comet_out_idxml = OPENMS_COMETADAPTER.out.idxml
    }

    // Run MS2PIP if specified
    if (params.use_ms2pip){
        MS2PIP(ch_comet_out_idxml.join(ch_clean_mzml_file))
        ch_versions = ch_versions.mix(MS2PIP.out.versions.ifEmpty(null))
        ch_comet_out_idxml_proceeding = MS2PIP.out.idxml
    } else {
        ch_comet_out_idxml_proceeding = ch_comet_out_idxml
    }

    // Write this information to an tsv file
    OPENMS_TEXTEXPORTER_COMET(ch_comet_out_idxml_proceeding)
    ch_versions = ch_versions.mix(OPENMS_COMETADAPTER.out.versions.ifEmpty(null))
    // Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(ch_comet_out_idxml_proceeding.join(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_PEPTIDEINDEXER.out.versions.ifEmpty(null))
    // Save indexed runs for later use to keep meta-run information. Sort based on file id
    ch_proceeding_idx = OPENMS_PEPTIDEINDEXER.out.idxml
            .map { meta, idxml ->
                [[id:meta.sample + "_" + meta.condition], meta, [id:meta.id, file:idxml]]
            }
            .groupTuple()
            .map {
                merge_id, meta, idxml ->
                // Create a custom comparator to sort based on id
                def comparator = { a, b -> a.id <=> b.id }
                [merge_id, meta.sort(comparator), idxml.sort(comparator).file]
            }

    ch_runs_to_merge = ch_proceeding_idx
            .map {
                merge_id, meta, idxml ->
                [merge_id, idxml]
            }
    // Merge aligned idXMLfiles
    OPENMS_IDMERGER(ch_runs_to_merge)
    ch_versions = ch_versions.mix(OPENMS_IDMERGER.out.versions.ifEmpty(null))

    // Extract PSM features for Percolator
    OPENMS_PSMFEATUREEXTRACTOR(OPENMS_IDMERGER.out.idxml)
    ch_versions = ch_versions.mix(OPENMS_PSMFEATUREEXTRACTOR.out.versions.ifEmpty(null))
    // Run Percolator
    OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
    ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions.ifEmpty(null))
    // Filter by percolator q-value
    OPENMS_IDFILTER_Q_VALUE(OPENMS_PERCOLATORADAPTER.out.idxml.flatMap { it -> [tuple(it[0], it[1], null)] })
    ch_versions = ch_versions.mix(OPENMS_IDFILTER_Q_VALUE.out.versions.ifEmpty(null))

    //
    // SUBWORKFLOW: Refine the FDR values on the predicted subset
    //
    if (params.refine_fdr_on_predicted_subset && params.predict_class_1) {
        // Run the following subworkflow
        REFINE_FDR (
            OPENMS_IDFILTER_Q_VALUE.out.idxml,
            OPENMS_PSMFEATUREEXTRACTOR.out.idxml,
            peptides_class_1_alleles
        )
        ch_versions = ch_versions.mix(REFINE_FDR.out.versions.ifEmpty(null))
        // Define the outcome of the paramer to a fixed variable
        filter_q_value = REFINE_FDR.out.filter_refined_q_value
    } else {
        // Make sure that the columns that consists of the ID's, sample names and the idXML file names are returned
        filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.idxml
    }

    //
    // SUBWORKFLOW: QUANT
    //
    if (!params.skip_quantification) {
        // Rip post-percolator idXML files
        OPENMS_IDRIPPER( OPENMS_PERCOLATORADAPTER.out.idxml ).ripped
                .join(filter_q_value)
                .join(ch_proceeding_idx)
                // [merge_id, [idxml_run1, idxml_run2, ...], pout_filtered, meta, [pre_perc_run1_idxml, pre_perc_run2_idxml, ...]]
                .map {merge_id, ripped, fdrfiltered, meta, indexed -> [meta, ripped, fdrfiltered] }
                .transpose()
                // [meta_run1, idxml_run1, pout_filtered] [meta_run2, idxml_run2, pout_filtered] ...
                .set { ch_runs_to_be_filtered }

        // Filter runs based on fdr filtered coprocessed percolator output.
        // NOTE: This is an alternative filtering method that will be replaced by IDFilter with new release of OpenMS
        PYOPENMS_IDFILTER( ch_runs_to_be_filtered ).filtered
                .map { meta, idxml ->
                        [[id:meta.sample + "_" + meta.condition], meta, [id:meta.id, file:idxml]]
                    }
                .groupTuple()
                .map { merge_id, meta, idxml ->
                    // Create a custom comparator to sort based on id
                    def comparator = { a, b -> a.id <=> b.id }
                    [merge_id, meta.sort(comparator), idxml.sort(comparator).file]
                    }
                .set { ch_runs_to_be_aligned }

        MAP_ALIGNMENT(
            ch_runs_to_be_aligned,
            ch_clean_mzml_file
        )
        ch_versions = ch_versions.mix(MAP_ALIGNMENT.out.versions.ifEmpty(null))

        // We need to merge groupwise the aligned idxml files together to use them as id_ext in featurefinder
        OPENMS_IDMERGER_QUANT(MAP_ALIGNMENT.out.aligned_idxml
                                    .map { meta, aligned_idxml ->
                                        [[id:meta.sample + "_" + meta.condition], aligned_idxml] }
                                    .groupTuple())

        // Manipulate channels such that we end up with : [meta, mzml, run_idxml, merged_runs_idxml]
        MAP_ALIGNMENT.out.aligned_mzml
                .join(MAP_ALIGNMENT.out.aligned_idxml)
                .map {meta, mzml, idxml ->
                    [[id:meta.sample + "_" + meta.condition], meta, [id:meta.id, file:mzml], [id:meta.id, file:idxml]] }
                .groupTuple()
                .map { merge_id, meta, mzml, idxml ->
                    // Create a custom comparator to sort based on id
                    def comparator = { a, b -> a.id <=> b.id }
                    [merge_id, meta.sort(comparator), mzml.sort(comparator).file, idxml.sort(comparator).file]}
                .join(OPENMS_IDMERGER_QUANT.out.idxml)
                .map { merge_id, meta, mzml, idxml, merged_idxml ->
                    [meta, mzml, idxml, merged_idxml] }
                .transpose()
                .set { ch_runs_to_be_quantified }

        PROCESS_FEATURE (
            ch_runs_to_be_quantified
        )
        ch_versions = ch_versions.mix(PROCESS_FEATURE.out.versions.ifEmpty(null))

    } else {
        // Prepare for check if file is empty
        OPENMS_TEXTEXPORTER_ID(OPENMS_IDFILTER_Q_VALUE.out.idxml)
        // Return an error message when there is only a header present in the document
        OPENMS_TEXTEXPORTER_ID.out.tsv.map {
            meta, tsv -> if (tsv.size() < 130) {
            log.warn "It seems that there were no significant hits found for this sample: " + meta.sample + "\nPlease consider incrementing the '--fdr_threshold' after removing the work directory or to exclude this sample. "
            }
        }
        OPENMS_MZTABEXPORTER_ID(filter_q_value)
    }

    //
    // SUBWORKFLOW: Predict class I (neoepitopes)
    //
    if (params.predict_class_1 & !params.skip_quantification) {
        PREDICT_CLASS1 (
            PROCESS_FEATURE.out.mztab,
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
    if (params.predict_class_2 & !params.skip_quantification) {
        PREDICT_CLASS2 (
            PROCESS_FEATURE.out.mztab,
            peptides_class_2_alleles,
            ch_vcf_from_sheet
        )
        ch_versions = ch_versions.mix(PREDICT_CLASS2.out.versions.ifEmpty(null))
        ch_predicted_possible_neoepitopes_II = PREDICT_CLASS2.out.ch_predicted_possible_neoepitopes
    } else {
        ch_predicted_possible_neoepitopes_II = Channel.empty()
    }

    if (params.annotate_ions) {
        // Alter the annotation of the filtered q value
        ch_filtered_idxml = filter_q_value.map { ident, meta, idxml -> [meta.id, idxml] }
        // Join the ch_filtered_idxml and the ch_mzml_file
        ch_raw_spectra_data = ch_clean_mzml_file.map {meta, mzml -> [meta.sample + '_' + meta.condition, mzml] }
            .groupTuple()
            .join(ch_filtered_idxml)
        // Annotate spectra with ion fragmentation information
        PYOPENMS_IONANNOTATOR(ch_raw_spectra_data)
        ch_versions = ch_versions.mix(PYOPENMS_IONANNOTATOR.out.versions.ifEmpty(null))
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
        workflow_summary = WorkflowMhcquant.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMhcquant.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.collect().ifEmpty([]),
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([])
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
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
