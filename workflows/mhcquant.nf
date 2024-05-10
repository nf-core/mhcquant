/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { TDF2MZML                   } from '../modules/local/tdf2mzml'
include { OPENMS_FILEFILTER          } from '../modules/local/openms_filefilter'
include { OPENMS_COMETADAPTER        } from '../modules/local/openms_cometadapter'
include { OPENMS_PEPTIDEINDEXER      } from '../modules/local/openms_peptideindexer'
include { MS2RESCORE                 } from '../modules/local/ms2rescore'
include { OPENMS_PSMFEATUREEXTRACTOR } from '../modules/local/openms_psmfeatureextractor'
include { OPENMS_PERCOLATORADAPTER   } from '../modules/local/openms_percolatoradapter'
include { PYOPENMS_IONANNOTATOR      } from '../modules/local/pyopenms_ionannotator'
include { OPENMS_TEXTEXPORTER        } from '../modules/local/openms_textexporter'
include { OPENMS_MZTABEXPORTER       } from '../modules/local/openms_mztabexporter'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { INCLUDE_PROTEINS } from '../subworkflows/local/include_proteins'
include { REFINE_FDR       } from '../subworkflows/local/refine_fdr'
include { QUANT            } from '../subworkflows/local/quant'
include { PREDICT_CLASS1   } from '../subworkflows/local/predict_class1'
include { PREDICT_CLASS2   } from '../subworkflows/local/predict_class2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { THERMORAWFILEPARSER                        } from '../modules/nf-core/thermorawfileparser/main'
include { OPENMS_DECOYDATABASE                       } from '../modules/nf-core/openms/decoydatabase/main'
include { OPENMS_PEAKPICKERHIRES                     } from '../modules/nf-core/openms/peakpickerhires/main'
include { OPENMS_IDMERGER                            } from '../modules/nf-core/openms/idmerger/main'
include { OPENMS_IDSCORESWITCHER                     } from '../modules/nf-core/openms/idscoreswitcher/main.nf'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE } from '../modules/nf-core/openms/idfilter/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                           } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_mhcquant_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MHCQUANT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta       // channel: reference database read in from --fasta

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
        .branch {
            meta, filename ->
                raw : meta.ext == 'raw'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                mzml : meta.ext == 'mzml'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                tdf : meta.ext == 'd'
                    return [ meta.subMap('id', 'sample', 'condition', 'group_count', 'spectra'), filename ]
                other : true }
        .set { branched_ms_files }

    //
    // SUBWORKFLOW: Include protein information
    //
    // TODO: Temporary disabled because of outdated vcf parsing
    //if (params.include_proteins_from_vcf) {
    //    // Include the proteins from the vcf file to the fasta file
    //    INCLUDE_PROTEINS(ch_fasta)
    //    ch_versions = ch_versions.mix(INCLUDE_PROTEINS.out.versions)
    //    ch_fasta_file = INCLUDE_PROTEINS.out.ch_fasta_file
    //    ch_vcf_from_sheet = INCLUDE_PROTEINS.out.ch_vcf_from_sheet
    //} else {
    //    ch_fasta_file = ch_fasta
    //    ch_vcf_from_sheet = Channel.empty()
    //}
    if (!params.skip_decoy_generation) {
        // Generate reversed decoy database
        OPENMS_DECOYDATABASE(ch_fasta)
        ch_versions = ch_versions.mix(OPENMS_DECOYDATABASE.out.versions)
        ch_decoy_db = OPENMS_DECOYDATABASE.out.decoy_fasta
                                .map{ meta, fasta -> [fasta] }
    } else {
        ch_decoy_db = ch_fasta.map{ meta, fasta -> [fasta] }
    }

    // Initialize channel for ms files that do not need to be converted
    ch_ms_files = branched_ms_files.mzml

    // Raw file conversion
    THERMORAWFILEPARSER(branched_ms_files.raw)
    ch_versions = ch_versions.mix(THERMORAWFILEPARSER.out.versions)
    ch_ms_files = ch_ms_files.mix(THERMORAWFILEPARSER.out.spectra)

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
    //    OPENMS_COMETADAPTER(ch_clean_mzml_file.combine(ch_fasta.map{ meta, fasta -> [fasta] }))
    //}
    OPENMS_COMETADAPTER(ch_clean_mzml_file.combine(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_COMETADAPTER.out.versions)

    // Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(OPENMS_COMETADAPTER.out.idxml.combine(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_PEPTIDEINDEXER.out.versions)

    // Save indexed runs for later use to keep meta-run information. Sort based on file id
    OPENMS_PEPTIDEINDEXER.out.idxml
        .map { meta, idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), meta] }
        .groupTuple()
        .set { merge_meta_map }

    OPENMS_PEPTIDEINDEXER.out.idxml
        .map { meta, idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), idxml] }
        .groupTuple()
        .set { ch_runs_to_merge }

    // Merge aligned idXMLfiles
    OPENMS_IDMERGER(ch_runs_to_merge)
    ch_versions = ch_versions.mix(OPENMS_IDMERGER.out.versions)

    // Run MS2Rescore
    ch_clean_mzml_file
            .map { meta, mzml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), mzml] }
            .groupTuple()
            .join(OPENMS_IDMERGER.out.idxml)
            .map { meta, mzml, idxml -> [meta, idxml, mzml, []] }
            .set { ch_ms2rescore_in }

    MS2RESCORE(ch_ms2rescore_in)
    ch_versions = ch_versions.mix(MS2RESCORE.out.versions)

    if (params.rescoring_engine == 'percolator') {
        // Extract PSM features for Percolator
        OPENMS_PSMFEATUREEXTRACTOR(MS2RESCORE.out.idxml.join(MS2RESCORE.out.feature_names))
        ch_versions = ch_versions.mix(OPENMS_PSMFEATUREEXTRACTOR.out.versions)

        // Run Percolator
        OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions)
        ch_rescored_runs = OPENMS_PERCOLATORADAPTER.out.idxml
    } else {
        log.warn "The rescoring engine is set to mokapot. This rescoring engine currently only supports psm-level-fdr via ms2rescore."
        // Switch comet e-value to mokapot q-value
        OPENMS_IDSCORESWITCHER(MS2RESCORE.out.idxml)
        ch_versions = ch_versions.mix(OPENMS_IDSCORESWITCHER.out.versions)
        ch_rescored_runs = OPENMS_IDSCORESWITCHER.out.idxml
    }

    // Filter by percolator q-value
    OPENMS_IDFILTER_Q_VALUE(ch_rescored_runs.map {group_meta, idxml -> [group_meta, idxml, []]})
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
        ch_filter_q_value = REFINE_FDR.out.filter_refined_q_value
    } else {
        // Make sure that the columns that consists of the ID's, sample names and the idXML file names are returned
        ch_filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.filtered
    }

    //
    // SUBWORKFLOW: QUANT
    //
    if (!params.skip_quantification) {
        QUANT(merge_meta_map, ch_rescored_runs, ch_filter_q_value, ch_clean_mzml_file)
        ch_versions = ch_versions.mix(QUANT.out.versions)
        ch_output = QUANT.out.consensusxml
    } else {
        ch_output = ch_filter_q_value
    }

    if (params.annotate_ions) {
        // Join the ch_filtered_idxml and the ch_mzml_file
        ch_clean_mzml_file.map { meta, mzml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), mzml] }
            .groupTuple()
            .join(ch_filter_q_value)
            .set{ ch_ion_annotator_input }

        // Annotate spectra with ion fragmentation information
        PYOPENMS_IONANNOTATOR( ch_ion_annotator_input )
        ch_versions = ch_versions.mix(PYOPENMS_IONANNOTATOR.out.versions)
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

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
