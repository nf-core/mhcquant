/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { PYOPENMS_CHROMATOGRAMEXTRACTOR } from '../modules/local/pyopenms/chromatogramextractor'
include { PYOPENMS_IONANNOTATOR          } from '../modules/local/pyopenms/ionannotator'
include { OPENMS_IDMASSACCURACY          } from '../modules/local/openms/idmassaccuracy/main'
include { OPENMS_TEXTEXPORTER            } from '../modules/nf-core/openms/textexporter/main'
include { SUMMARIZE_RESULTS              } from '../modules/local/pyopenms/summarize_results'
include { EPICORE                        } from '../modules/local/epicore'

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
    ch_sdrf        // channel: raw SDRF file, empty for plain samplesheet input
    ch_accession   //    val: resolved PRIDE accession, '' for plain samplesheet input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_multiqc_files = channel.empty()

    // Prepare spectra files (Decompress archives, convert to mzML, centroid if specified)
    PREPARE_SPECTRA(ch_samplesheet)

    // Decoy Database creation
    if (!params.skip_decoy_generation) {
        // Generate reversed decoy database
        OPENMS_DECOYDATABASE(ch_fasta)
        ch_decoy_db = OPENMS_DECOYDATABASE.out.decoy_fasta
    } else {
        ch_decoy_db = ch_fasta
    }

    // Optionally clean up mzML files
    if (params.filter_mzml){
        OPENMS_FILEFILTER(PREPARE_SPECTRA.out.mzml)
        ch_clean_mzml_file = OPENMS_FILEFILTER.out.mzml
    } else {
        ch_clean_mzml_file = PREPARE_SPECTRA.out.mzml
    }

    // Compute MS1 TICs for QC
    PYOPENMS_CHROMATOGRAMEXTRACTOR(ch_clean_mzml_file)
    ch_multiqc_files = ch_multiqc_files.mix(PYOPENMS_CHROMATOGRAMEXTRACTOR.out.csv.map{ meta, mzml -> mzml })

    // Prepare the comet input channel with global fasta or per-sample_condition fasta
    ch_comet_in = params.fasta ?
        ch_clean_mzml_file.combine(ch_decoy_db.map{ meta, fasta -> [fasta] }) :
        ch_clean_mzml_file
            .map { meta, mzml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), meta, mzml] }
            .combine(ch_decoy_db, by: 0)
            .map { groupKey, meta, mzml, fasta -> [meta, mzml, fasta] }

    // Run comet database search and index decoy and target hits
    OPENMSTHIRDPARTY_COMETADAPTER(ch_comet_in)

    // Prepare the peptideindexer channel with global fasta or per-sample_condition fasta
    ch_peptideindexer_in = params.fasta ?
        OPENMSTHIRDPARTY_COMETADAPTER.out.idxml.combine(ch_decoy_db.map{ meta, fasta -> [fasta] }) :
        OPENMSTHIRDPARTY_COMETADAPTER.out.idxml
            .map { meta, idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), meta, idxml] }
            .combine(ch_decoy_db, by: 0)
            .map { groupKey, meta, idxml, fasta -> [meta, idxml, fasta] }

    OPENMS_PEPTIDEINDEXER(ch_peptideindexer_in)

    // Compute mass errors for multiQC report
    OPENMS_IDMASSACCURACY(PREPARE_SPECTRA.out.mzml.join(OPENMS_PEPTIDEINDEXER.out.indexed_idxml))
    ch_multiqc_files = ch_multiqc_files.mix(OPENMS_IDMASSACCURACY.out.frag_err.map{ meta, frag_err -> frag_err })

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

    // Run MS2Rescore
    ch_clean_mzml_file
            .map { meta, mzml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), meta, mzml] }
            .groupTuple()
            .join(OPENMS_IDMERGER.out.idxml)
            .map { gkey, metas, mzmls, idxml ->
                // All replicates in a group share the same search params
                def meta = metas[0]
                def searchMeta = meta.subMap(meta.keySet() - ['id', 'sample', 'condition', 'group_count', 'spectra'])
                [[id: "${meta.sample}_${meta.condition}"] + searchMeta, idxml, mzmls, []]
            }
            .set { ch_rescore_in }

    //
    // SUBWORKFLOW: RESCORE WITH MOKKAPOT OR PERCOLATOR AND FILTER BY Q-VALUE ON LOCAL/GLOBAL FDR
    //
    RESCORE( ch_rescore_in, ch_multiqc_files )
    ch_multiqc_files = ch_multiqc_files.mix(RESCORE.out.multiqc_files)

    RESCORE.out.fdr_filtered_empty.subscribe { item ->
        log.warn """\
            No peptides passed FDR filtering for sample '${item[0].id}'.
            An empty TSV will be written.
            Consider raising '--fdr_threshold' or excluding this sample.
            """.stripIndent()
    }

    // GENERATE SPECTRUM LIBRARY
    if (params.generate_speclib) {
        OPENMSTHIRDPARTY_COMETADAPTER.out.idxml
                .map { meta, idxml -> [ [id: "${meta.sample}_${meta.condition}"], meta, idxml] }
                .combine(RESCORE.out.fdr_filtered, by:0)
                .map { groupKey, meta, comet_idxml, fdr_filtered_idxml -> [meta, comet_idxml, fdr_filtered_idxml] }
                .set { ch_fdrfilter_comet_idxml }

        // Backfilter Comet identifications with FDR threshold
        OPENMS_IDFILTER_FOR_SPECLIB(ch_fdrfilter_comet_idxml)

        //
        // SUBWORKFLOW: SPECLIB
        //
        SPECLIB(OPENMS_IDFILTER_FOR_SPECLIB.out.filtered, ch_clean_mzml_file)
    }

    //
    // SUBWORKFLOW: QUANT
    //
    if (params.quantify) {
        QUANT(merge_meta_map, RESCORE.out.rescored_runs, RESCORE.out.fdr_filtered, ch_clean_mzml_file, ch_sdrf, ch_accession)
        ch_output = QUANT.out.consensusxml.mix(RESCORE.out.fdr_filtered_empty)
    } else {
        ch_output = RESCORE.out.fdr_filtered.mix(RESCORE.out.fdr_filtered_empty)
    }

    // Annotate Ions for follow-up spectrum validation
    if (params.annotate_ions) {
        // Join mzml files and fdr-filtered idxml, reconstructing search params from ch_clean_mzml_file
        ch_clean_mzml_file.map { meta, mzml ->
                [ [id: "${meta.sample}_${meta.condition}"],
                  meta.subMap(meta.keySet() - ['id', 'sample', 'condition', 'group_count', 'spectra']),
                  mzml ]
            }
            .groupTuple()
            .join(RESCORE.out.fdr_filtered)
            .map { key, searchMetas, mzmls, idxml -> [key + searchMetas[0], mzmls, idxml] }
            .set{ ch_ion_annotator_input }

        // Annotate spectra with ion fragmentation information
        PYOPENMS_IONANNOTATOR( ch_ion_annotator_input )
    }

    OPENMS_TEXTEXPORTER(ch_output)

    // Process the tsv file to facilitate visualization with MultiQC.
    // Under --quantify, attach each group's per-run trafoXMLs so alignment residuals can be plotted.
    if (params.quantify) {
        ch_summarize_input = OPENMS_TEXTEXPORTER.out.tsv
            .map { meta, tsv -> [meta.id, meta, tsv] }
            .join( QUANT.out.trafoxml.map { meta, trafoxml -> [meta.id, trafoxml] }, remainder: true )
            .map { _id, meta, tsv, trafoxml -> [meta, tsv, trafoxml ?: []] }
    } else {
        ch_summarize_input = OPENMS_TEXTEXPORTER.out.tsv.map { meta, tsv -> [meta, tsv, []] }
    }
    SUMMARIZE_RESULTS(ch_summarize_input)

    //
    // EPICORE
    //
    if (params.epicore) {
        EPICORE(ch_fasta.map{ fasta -> fasta.last()}.first(), SUMMARIZE_RESULTS.out.epicore_input)
        ch_multiqc_files = ch_multiqc_files.mix(
            EPICORE.out.length_dist,
            EPICORE.out.intensity_hist
        )
    }

    //
    // Collate MultiQC files
    //
    ch_multiqc_files = ch_multiqc_files.mix(
        SUMMARIZE_RESULTS.out.hist_mz,
        SUMMARIZE_RESULTS.out.hist_rt,
        SUMMARIZE_RESULTS.out.hist_scores,
        SUMMARIZE_RESULTS.out.hist_im,
        SUMMARIZE_RESULTS.out.xcorr,
        SUMMARIZE_RESULTS.out.lengths,
        SUMMARIZE_RESULTS.out.intensities,
        SUMMARIZE_RESULTS.out.rt_calibration,
        SUMMARIZE_RESULTS.out.aligned_residuals,
        params.epicore ? EPICORE.out.stats : SUMMARIZE_RESULTS.out.epicore_input.map { meta, tsv, stats -> stats }
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'mhcquant_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'mhcquant'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
