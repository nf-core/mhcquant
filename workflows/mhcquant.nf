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
if (params.fasta)   { params.fasta = params.fasta       }

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
include { THERMORAWFILEPARSER }                                             from '../modules/local/thermorawfileparser'
include { TDF2MZML }                                                        from '../modules/local/tdf2mzml'
include { OPENMS_FILEFILTER }                                               from '../modules/local/openms_filefilter'
include { OPENMS_COMETADAPTER }                                             from '../modules/local/openms_cometadapter'
include { OPENMS_PEPTIDEINDEXER }                                           from '../modules/local/openms_peptideindexer'
include { MS2RESCORE }                                                      from '../modules/local/ms2rescore'
include { OPENMS_IDSCORESWITCHER }                                          from '../modules/local/openms_idscoreswitcher'

include { OPENMS_PSMFEATUREEXTRACTOR }                                      from '../modules/local/openms_psmfeatureextractor'
include { OPENMS_PERCOLATORADAPTER }                                        from '../modules/local/openms_percolatoradapter'
include { PYOPENMS_IONANNOTATOR }                                           from '../modules/local/pyopenms_ionannotator'

include { OPENMS_TEXTEXPORTER }                                             from '../modules/local/openms_textexporter'
include { OPENMS_MZTABEXPORTER }                                            from '../modules/local/openms_mztabexporter'


//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { INCLUDE_PROTEINS }                                                from '../subworkflows/local/include_proteins'
include { REFINE_FDR }                                                      from '../subworkflows/local/refine_fdr'
include { QUANT }                                                           from '../subworkflows/local/quant'
include { PREDICT_CLASS1 }                                                  from '../subworkflows/local/predict_class1'
include { PREDICT_CLASS2 }                                                  from '../subworkflows/local/predict_class2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { OPENMS_DECOYDATABASE                       } from '../modules/nf-core/openms/decoydatabase/main'
include { OPENMS_PEAKPICKERHIRES                     } from '../modules/nf-core/openms/peakpickerhires/main'
include { OPENMS_IDMERGER                            } from '../modules/nf-core/openms/idmerger/main'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE } from '../modules/nf-core/openms/idfilter/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []
// Sort closure for merging and splitting files
def sortById = { a, b -> a.id <=> b.id }


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
                    return [ meta.subMap('id', 'sample', 'condition'), filename ]
                mzml : meta.ext == 'mzml'
                    return [ meta.subMap('id', 'sample', 'condition'), filename ]
                tdf : meta.ext == 'd'
                    return [ meta.subMap('id', 'sample', 'condition'), filename ]
                other : true }
        .set { branched_ms_files }

    // Input fasta file
    Channel.fromPath(params.fasta)
        .map{ fasta -> [[id:fasta.getBaseName()], fasta] }
        .ifEmpty { error ("params.fasta was empty - no input file supplied") }
        .set { fasta_file }

    //
    // SUBWORKFLOW: Include protein information
    //
    // TODO: Temporary disabled because of outdated vcf parsing
    //if (params.include_proteins_from_vcf) {
    //    // Include the proteins from the vcf file to the fasta file
    //    INCLUDE_PROTEINS(fasta_file)
    //    ch_versions = ch_versions.mix(INCLUDE_PROTEINS.out.versions)
    //    ch_fasta_file = INCLUDE_PROTEINS.out.ch_fasta_file
    //    ch_vcf_from_sheet = INCLUDE_PROTEINS.out.ch_vcf_from_sheet
    //} else {
    //    ch_fasta_file = fasta_file
    //    ch_vcf_from_sheet = Channel.empty()
    //}
    if (!params.skip_decoy_generation) {
        // Generate reversed decoy database
        OPENMS_DECOYDATABASE(fasta_file)
        ch_versions = ch_versions.mix(OPENMS_DECOYDATABASE.out.versions)
        ch_decoy_db = OPENMS_DECOYDATABASE.out.decoy_fasta
                                .map{ meta, fasta -> [fasta] }
    } else {
        ch_decoy_db = fasta_file.map{ meta, fasta -> [fasta] }
    }

    // If mzml files are specified, they are encapsulated in a list [meta, [mzml]]. We need to extract the path for grouping later
    ch_ms_files = branched_ms_files.mzml.map{ meta, mzml -> [meta, mzml[0]]}
    // Raw file conversion
    THERMORAWFILEPARSER(branched_ms_files.raw)
    ch_versions = ch_versions.mix(THERMORAWFILEPARSER.out.versions)
    ch_ms_files = ch_ms_files.mix(THERMORAWFILEPARSER.out.mzml)

    // timsTOF data conversion
    TDF2MZML(branched_ms_files.tdf)
    ch_versions = ch_versions.mix(TDF2MZML.out.versions)
    ch_ms_files = ch_ms_files.mix(TDF2MZML.out.mzml)

    // Optional: Run Peak Picking as Preprocessing
    if (params.run_centroidisation) {
        OPENMS_PEAKPICKERHIRES(ch_ms_files)
        ch_versions = ch_versions.mix(OPENMS_PEAKPICKERHIRES.out.versions)
        ch_mzml_file = OPENMS_PEAKPICKERHIRES.out.mzml
    } else {
        ch_mzml_file = ch_ms_files
    }

    // Optionally clean up mzML files
    if (params.filter_mzml){
        OPENMS_FILEFILTER(ch_mzml_file)
        ch_versions = ch_versions.mix(OPENMS_FILEFILTER.out.versions)
        ch_clean_mzml_file = OPENMS_FILEFILTER.out.cleaned_mzml
    } else {
        ch_clean_mzml_file = ch_mzml_file
    }

    // Run comet database search
    // TODO: Fix accordingly with vcf parsing
    //if (params.include_proteins_from_vcf) {
    //    OPENMS_COMETADAPTER(ch_clean_mzml_file.join(ch_decoy_db, remainder:true))
    //} else {
    //    OPENMS_COMETADAPTER(ch_clean_mzml_file.combine(ch_fasta_file.map{ meta, fasta -> [fasta] }))
    //}
    OPENMS_COMETADAPTER(ch_clean_mzml_file.combine(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_COMETADAPTER.out.versions)

    // Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(OPENMS_COMETADAPTER.out.idxml.combine(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_PEPTIDEINDEXER.out.versions)

    // Save indexed runs for later use to keep meta-run information. Sort based on file id
    OPENMS_PEPTIDEINDEXER.out.idxml
            .map { meta, idxml -> [[id: meta.sample + '_' + meta.condition], meta] }
            .groupTuple( sort: sortById )
            .set { merge_meta_map }

    OPENMS_PEPTIDEINDEXER.out.idxml
            .map { meta, idxml -> [[id: meta.sample + '_' + meta.condition], idxml] }
            .groupTuple()
            .set { ch_runs_to_merge }

    // Merge aligned idXMLfiles
    OPENMS_IDMERGER(ch_runs_to_merge)
    ch_versions = ch_versions.mix(OPENMS_IDMERGER.out.versions)

    // Run MS2Rescore
    ch_clean_mzml_file
            .map { meta, mzml -> [[id: meta.sample + '_' + meta.condition], mzml] }
            .groupTuple()
            .join(OPENMS_IDMERGER.out.idxml)
            .map { meta, mzml, idxml -> [meta, idxml, mzml, []] }
            .set { ch_ms2rescore_in }

    MS2RESCORE(ch_ms2rescore_in)
    ch_versions = ch_versions.mix(MS2RESCORE.out.versions)

    if (params.rescoring_engine == 'percolator') {
        // Extract PSM features for Percolator
        OPENMS_PSMFEATUREEXTRACTOR(MS2RESCORE.out.idxml
                                        .join(MS2RESCORE.out.feature_names))
        ch_versions = ch_versions.mix(OPENMS_PSMFEATUREEXTRACTOR.out.versions)

        // Run Percolator
        OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions)
        ch_rescored_runs = OPENMS_PERCOLATORADAPTER.out.idxml
    } else {
        log.warn "The rescoring engine is set to mokapot. This rescoring engine currently only supports psm-level-fdr via ms2rescore."
        // TODO: remove whitelist argument from idscoreswitcher
        OPENMS_IDSCORESWITCHER(MS2RESCORE.out.idxml
                                    .map { meta, idxml -> [meta, idxml, []] })
        ch_rescored_runs = OPENMS_IDSCORESWITCHER.out.switched_idxml.map { tuple -> tuple.findAll { it != [] }}
    }

    // Filter by percolator q-value
    OPENMS_IDFILTER_Q_VALUE(ch_rescored_runs.flatMap { it -> [tuple(it[0], it[1], [])] })
    ch_versions = ch_versions.mix(OPENMS_IDFILTER_Q_VALUE.out.versions)

    //
    // SUBWORKFLOW: Refine the FDR values on the predicted subset
    //
    if (params.refine_fdr_on_predicted_subset && params.predict_class_1) {
        // Run the following subworkflow
        REFINE_FDR (
            OPENMS_IDFILTER_Q_VALUE.out.filtered,
            OPENMS_PSMFEATUREEXTRACTOR.out.idxml,
            peptides_class_1_alleles
        )
        ch_versions = ch_versions.mix(REFINE_FDR.out.versions)
        // Define the outcome of the paramer to a fixed variable
        filter_q_value = REFINE_FDR.out.filter_refined_q_value
    } else {
        // Make sure that the columns that consists of the ID's, sample names and the idXML file names are returned
        filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.filtered
    }

    //
    // SUBWORKFLOW: QUANT
    //
    if (!params.skip_quantification) {
        QUANT(merge_meta_map, ch_rescored_runs, filter_q_value, ch_clean_mzml_file)
        ch_versions = ch_versions.mix(QUANT.out.versions)
        ch_output = QUANT.out.consensusxml
    } else {
        ch_output = filter_q_value
    }

    // Prepare for check if file is empty
    OPENMS_TEXTEXPORTER(ch_output)
    ch_versions = ch_versions.mix(OPENMS_TEXTEXPORTER.out.versions)
    // Return an error message when there is only a header present in the document
    OPENMS_TEXTEXPORTER.out.tsv.map {
        meta, tsv -> if (tsv.size() < 130) {
        log.warn "It seems that there were no significant hits found for this sample: " + meta.sample + "\nPlease consider incrementing the '--fdr_threshold' after removing the work directory or to exclude this sample. "
        }
    }

    OPENMS_MZTABEXPORTER(ch_output)
    ch_versions = ch_versions.mix(OPENMS_MZTABEXPORTER.out.versions)

    //
    // SUBWORKFLOW: Predict class I (neoepitopes)
    //
    // TODO: Temporary disabled because of outdated vcf parsing
    //if (params.predict_class_1 & !params.skip_quantification) {
    //    PREDICT_CLASS1 (
    //        OPENMS_MZTABEXPORTER.out.mztab,
    //        peptides_class_1_alleles,
    //        ch_vcf_from_sheet
    //    )
    //    ch_versions = ch_versions.mix(PREDICT_CLASS1.out.versions)
    //    ch_predicted_possible_neoepitopes = PREDICT_CLASS1.out.ch_predicted_possible_neoepitopes
    //} else {
    //    ch_predicted_possible_neoepitopes = Channel.empty()
    //}
    //
    ////
    //// SUBWORKFLOW: Predict class II (neoepitopes)
    ////
    //if (params.predict_class_2 & !params.skip_quantification) {
    //    PREDICT_CLASS2 (
    //        OPENMS_MZTABEXPORTER.out.mztab,
    //        peptides_class_2_alleles,
    //        ch_vcf_from_sheet
    //    )
    //    ch_versions = ch_versions.mix(PREDICT_CLASS2.out.versions)
    //    ch_predicted_possible_neoepitopes_II = PREDICT_CLASS2.out.ch_predicted_possible_neoepitopes
    //} else {
    //    ch_predicted_possible_neoepitopes_II = Channel.empty()
    //}

    if (params.annotate_ions) {
        // Join the ch_filtered_idxml and the ch_mzml_file
        ch_clean_mzml_file.map { meta, mzml -> [[id: meta.sample + '_' + meta.condition], mzml] }
            .groupTuple()
            .join(filter_q_value)
            .set{ ch_ion_annotator_input }

        // Annotate spectra with ion fragmentation information
        PYOPENMS_IONANNOTATOR( ch_ion_annotator_input )
        ch_versions = ch_versions.mix(PYOPENMS_IONANNOTATOR.out.versions)
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
    NfcoreTemplate.dump_parameters(workflow, params)
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
