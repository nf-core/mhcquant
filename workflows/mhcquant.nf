/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { PYOPENMS_CHROMATOGRAMEXTRACTOR } from '../modules/local/pyopenms/chromatogramextractor'
include { DATAMASH_HISTOGRAM             } from '../modules/local/datamash_histogram'
include { PYOPENMS_IONANNOTATOR          } from '../modules/local/pyopenms/ionannotator'
include { OPENMS_TEXTEXPORTER            } from '../modules/local/openms/textexporter'
include { SUMMARIZE_RESULTS              } from '../modules/local/pyopenms/summarize_results'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { PREPARE_SPECTRA } from '../subworkflows/local/prepare_spectra'
include { RESCORE         } from '../subworkflows/local/rescore'
include { SPECLIB         } from '../subworkflows/local/speclib'
include { QUANT           } from '../subworkflows/local/quant'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { OPENMS_FILEFILTER                              } from '../modules/nf-core/openms/filefilter/main'
include { OPENMS_DECOYDATABASE                           } from '../modules/nf-core/openms/decoydatabase/main'
include { OPENMS_IDMASSACCURACY                          } from '../modules/nf-core/openms/idmassaccuracy/main'
include { OPENMSTHIRDPARTY_COMETADAPTER                  } from '../modules/nf-core/openmsthirdparty/cometadapter/main'
include { OPENMS_PEPTIDEINDEXER                          } from '../modules/nf-core/openms/peptideindexer/main'
include { OPENMS_IDMERGER                                } from '../modules/nf-core/openms/idmerger/main'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_FOR_SPECLIB } from '../modules/nf-core/openms/idfilter/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                         } from '../subworkflows/local/utils_nfcore_mhcquant_pipeline'

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

    // Prepare spectra files (Decompress archives, convert to mzML, centroid if specified)
    PREPARE_SPECTRA(ch_samplesheet)
    ch_versions = ch_versions.mix(PREPARE_SPECTRA.out.versions)

    // Decoy Database creation
    if (!params.skip_decoy_generation) {
        // Generate reversed decoy database
        OPENMS_DECOYDATABASE(ch_fasta)
        ch_versions = ch_versions.mix(OPENMS_DECOYDATABASE.out.versions)
        ch_decoy_db = OPENMS_DECOYDATABASE.out.decoy_fasta.map{ meta, fasta -> [fasta] }
    } else {
        ch_decoy_db = ch_fasta.map{ meta, fasta -> [fasta] }
    }

    // Optionally clean up mzML files
    if (params.filter_mzml){
        OPENMS_FILEFILTER(PREPARE_SPECTRA.out.mzml)
        ch_versions = ch_versions.mix(OPENMS_FILEFILTER.out.versions)
        ch_clean_mzml_file = OPENMS_FILEFILTER.out.mzml
    } else {
        ch_clean_mzml_file = PREPARE_SPECTRA.out.mzml
    }

    // Compute MS1 TICs for QC
    PYOPENMS_CHROMATOGRAMEXTRACTOR(ch_clean_mzml_file)
    ch_versions = ch_versions.mix(PYOPENMS_CHROMATOGRAMEXTRACTOR.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PYOPENMS_CHROMATOGRAMEXTRACTOR.out.csv.map{ meta, mzml -> mzml })

    // Run comet database search
    OPENMSTHIRDPARTY_COMETADAPTER(ch_clean_mzml_file.combine(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMSTHIRDPARTY_COMETADAPTER.out.versions)

    // Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(OPENMSTHIRDPARTY_COMETADAPTER.out.idxml.combine(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_PEPTIDEINDEXER.out.versions)

    // Compute mass errors for multiQC report
    OPENMS_IDMASSACCURACY(PREPARE_SPECTRA.out.mzml.join(OPENMS_PEPTIDEINDEXER.out.indexed_idxml))
    ch_versions = ch_versions.mix(OPENMS_IDMASSACCURACY.out.versions)
    // Bin and count mass errors for multiQC report
    DATAMASH_HISTOGRAM(OPENMS_IDMASSACCURACY.out.frag_err)
    ch_versions = ch_versions.mix(DATAMASH_HISTOGRAM.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DATAMASH_HISTOGRAM.out.binned_tsv.map{ meta, frag_err_hist -> frag_err_hist })

    // Save indexed runs for later use to keep meta-run information. Sort based on file id
    OPENMS_PEPTIDEINDEXER.out.indexed_idxml
        .map { meta, idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), meta] }
        .groupTuple()
        .set { merge_meta_map }

    OPENMS_PEPTIDEINDEXER.out.indexed_idxml
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
            .set { ch_rescore_in }

    //
    // SUBWORKFLOW: RESCORE WITH MOKKAPOT OR PERCOLATOR AND FILTER BY Q-VALUE ON LOCAL/GLOBAL FDR
    //
    RESCORE( ch_rescore_in, ch_multiqc_files )
    ch_versions = ch_versions.mix(RESCORE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(RESCORE.out.multiqc_files)

    // GENERATE SPECTRUM LIBRARY
    if (params.generate_speclib) {
        OPENMSTHIRDPARTY_COMETADAPTER.out.idxml
                .map { meta, idxml -> [ [id: "${meta.sample}_${meta.condition}"], meta, idxml] }
                .combine(RESCORE.out.fdr_filtered, by:0)
                .map { groupKey, meta, comet_idxml, fdr_filtered_idxml -> [meta, comet_idxml, fdr_filtered_idxml] }
                .set { ch_fdrfilter_comet_idxml }

        // Backfilter Comet identifications with FDR threshold
        OPENMS_IDFILTER_FOR_SPECLIB(ch_fdrfilter_comet_idxml)
            .filtered
            .filter { meta, idxml -> idxml.countLines() > 130 }
            .set { ch_fdrfilter_comet_idxml_filtered }

        //
        // SUBWORKFLOW: SPECLIB
        //
        SPECLIB(ch_fdrfilter_comet_idxml_filtered, ch_clean_mzml_file)
        ch_versions = ch_versions.mix(SPECLIB.out.versions)
    }

    //
    // SUBWORKFLOW: QUANT
    //
    if (params.quantify) {
        QUANT(merge_meta_map, RESCORE.out.rescored_runs, RESCORE.out.fdr_filtered, ch_clean_mzml_file)
        ch_versions = ch_versions.mix(QUANT.out.versions)
        ch_output = QUANT.out.consensusxml
    } else {
        ch_output = RESCORE.out.fdr_filtered
    }

    // Annotate Ions for follow-up spectrum validation
    if (params.annotate_ions) {
        // Join the ch_filtered_idxml and the ch_mzml_file
        ch_clean_mzml_file.map { meta, mzml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), mzml] }
            .groupTuple()
            .join(RESCORE.out.fdr_filtered)
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

    // Process the tsv file to facilitate visualization with MultiQC
    SUMMARIZE_RESULTS(OPENMS_TEXTEXPORTER.out.tsv)
    ch_versions = ch_versions.mix(SUMMARIZE_RESULTS.out.versions)

    ch_multiqc_files = ch_multiqc_files.mix(
        SUMMARIZE_RESULTS.out.hist_mz,
        SUMMARIZE_RESULTS.out.hist_rt,
        SUMMARIZE_RESULTS.out.hist_scores,
        SUMMARIZE_RESULTS.out.hist_xcorr,
        SUMMARIZE_RESULTS.out.lengths,
        SUMMARIZE_RESULTS.out.stats
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'mhcquant_software_'  + 'mqc_'  + 'versions.yml',
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
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

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
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
