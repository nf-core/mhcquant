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

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
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

// Variant
if (params.include_proteins_from_vcf)  {
    vcf_sheet = file(params.vcf_sheet, checkIfExists: true)
    Channel.from( vcf_sheet )
        .splitCsv(header: ['Sample', 'VCF_FileName'], sep:'\t', skip: 1)
        .map { col -> tuple("${col.Sample}", file("${col.VCF_FileName}"),) }
        .set { ch_vcf_from_sheet }
    }

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
fdr_adj_threshold = (params.fdr_threshold == '0.01') ? '0.05' : params.fdr_threshold

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def openms_map_aligner_identification_options = modules['openms_map_aligner_identification']
def openms_comet_adapter_options = modules['openms_comet_adapter']
def generate_proteins_from_vcf_options = modules['generate_proteins_from_vcf']
def percolator_adapter_options = modules['percolator_adapter']
def id_filter_options = modules['id_filter']
def id_filter_for_alignment_options = id_filter_options.clone()
def id_filter_whitelist_options = modules['id_filter_whitelist']

id_filter_options.args += " -score:pep " + params.fdr_threshold
id_filter_for_alignment_options.args += " -score:pep "  + fdr_adj_threshold
openms_comet_adapter_options.args += x_ions + z_ions + c_ions + a_ions + NL_ions + rm_precursor
generate_proteins_from_vcf_options.args += variant_indel_filter + variant_snp_filter + variant_frameshift_filter
percolator_adapter_options.args += fdr_level
percolator_adapter_options.suffix = "all_ids_merged_psm_perc"

def percolator_adapter_klammer_options = percolator_adapter_options.clone()
percolator_adapter_klammer_options.args += " -klammer"

def id_filter_qvalue_options = id_filter_options.clone()
id_filter_qvalue_options.suffix = "filtered"

include { hasExtension }                                    from '../modules/local/functions'

include { INPUT_CHECK }                                     from '../subworkflow/local/input_check'                          addParams( options: [:] )
include { GENERATE_PROTEINS_FROM_VCF }                      from '../modules/local/generate_proteins_from_vcf'                       addParams( options: generate_proteins_from_vcf_options )
include { OPENMS_DECOYDATABASE }                            from '../modules/local/openms_decoydatabase'                             addParams( options: [:] )
include { OPENMS_THERMORAWFILEPARSER }                      from '../modules/local/openms_thermorawfileparser'                       addParams( options: [:] )
include { OPENMS_PEAKPICKERHIRES }                          from '../modules/local/openms_peakpickerhires'                           addParams( options: [:] )
include { OPENMS_COMETADAPTER }                             from '../modules/local/openms_cometadapter'                              addParams( options: openms_comet_adapter_options )
include { OPENMS_PEPTIDEINDEXER }                           from '../modules/local/openms_peptideindexer'                            addParams( options: [:] )
include { OPENMS_FALSEDISCOVERYRATE }                       from '../modules/local/openms_falsediscoveryrate'                        addParams( options: [:] )
include { OPENMS_IDFILTER as OPENMS_IDFILTER_FOR_ALIGNMENT }from '../modules/local/openms_idfilter'                                  addParams( options: id_filter_for_alignment_options )
include { OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE }      from '../modules/local/openms_idfilter'                                  addParams( options: id_filter_qvalue_options )
include { OPENMS_MAPALIGNERIDENTIFICATION }                 from '../modules/local/openms_mapaligneridentification'                  addParams( options: openms_map_aligner_identification_options )

include {
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERMZML
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERIDXML }        from '../modules/local/openms_maprttransformer'                 addParams( options: [:] )

include { OPENMS_IDMERGER }                                 from '../modules/local/openms_idmerger'                                  addParams( options: [:] )
include { OPENMS_PSMFEATUREEXTRACTOR }                      from '../modules/local/openms_psmfeatureextractor'                       addParams( options: [:] )
include { OPENMS_PERCOLATORADAPTER }                        from '../modules/local/openms_percolatoradapter'                         addParams( options: percolator_adapter_options )
include { OPENMS_PERCOLATORADAPTER as OPENMS_PERCOLATORADAPTER_KLAMMER } from '../modules/local/openms_percolatoradapter'            addParams( options: percolator_adapter_klammer_options )

include { REFINE_FDR_ON_PREDICTED_SUBSET }                  from '../subworkflow/local/refine_fdr_on_predicted_subset'       addParams( run_percolator_options : percolator_adapter_options, filter_options: id_filter_options, whitelist_filter_options: id_filter_whitelist_options)

include { OPENMS_FEATUREFINDERIDENTIFICATION }              from '../modules/local/openms_featurefinderidentification'               addParams( options: [:] )
include { OPENMS_FEATURELINKERUNLABELEDKD }                 from '../modules/local/openms_featurelinkerunlabeledkd'                  addParams( options: [:] )
include { OPENMS_IDCONFLICTRESOLVER }                       from '../modules/local/openms_idconflictresolver'                        addParams( options: [:] )
include { OPENMS_TEXTEXPORTER }                             from '../modules/local/openms_textexporter'                              addParams( options: [:] )
include { OPENMS_MZTABEXPORTER }                            from '../modules/local/openms_mztabexporter'                             addParams( options: [:] )

include { PREDICT_PEPTIDES_MHCFLURRY_CLASS_1 }              from '../modules/local/predict_peptides_mhcflurry_class_1'               addParams( options: [:] )
include { PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 }          from '../modules/local/preprocess_peptides_mhcnuggets_class_2'           addParams( options: [:] )
include { PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2 }             from '../modules/local/predict_peptides_mhcnuggets_class_2'              addParams( options: [:] )
include { POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 }         from '../modules/local/postprocess_peptides_mhcnuggets_class_2'          addParams( options: [:] )
include { PREDICT_POSSIBLE_NEOEPITOPES }                    from '../modules/local/predict_possible_neoepitopes'                     addParams( options: [:] )
include { PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES }            from '../modules/local/predict_possible_class_2_neoepitopes'             addParams( options: [:] )
include { RESOLVE_FOUND_NEOEPITOPES }                       from '../modules/local/resolve_found_neoepitopes'                        addParams( options: [:] )
include { RESOLVE_FOUND_CLASS_2_NEOEPITOPES }               from '../modules/local/resolve_found_class_2_neoepitopes'                addParams( options: [:] )
include { PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1 }           from '../modules/local/predict_neoepitopes_mhcflurry_class_1'            addParams( options: [:] )
include { PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 }       from '../modules/local/preprocess_neoepitopes_mhcnuggets_class_2'        addParams( options: [:] )
include { PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2 }          from '../modules/local/predict_neoepitopes_mhcnuggets_class_2'           addParams( options: [:] )
include { POSTPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 }      from '../modules/local/postprocess_neoepitopes_mhcnuggets_class_2'       addParams( options: [:] )

include { OPENMS_RTMODEL }                                  from '../modules/local/openms_rtmodel'                                   addParams( options: [:] )
include { OPENMS_RTPREDICT as OPENMS_RTPREDICT_FOUND_PEPTIDES}      from '../modules/local/openms_rtpredict'                         addParams( options: [suffix:"_id_files_for_rt_prediction_RTpredicted"] )
include { OPENMS_RTPREDICT as OPENMS_RTPREDICT_NEOEPITOPES}         from '../modules/local/openms_rtpredict'                         addParams( options: [suffix:"_txt_file_for_rt_prediction_RTpredicted"] )

include { CUSTOM_DUMPSOFTWAREVERSIONS }                     from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'       addParams( options: [publish_files : ['_versions.yml':'']] )

////////////////////////////////////////////////////
/* --              CREATE CHANNELS             -- */
////////////////////////////////////////////////////
// params.summary_params = [:]
include { SAMPLESHEET_CHECK } from '../modules/local/samplesheet_check' addParams( )
////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow MHCQUANT {

    ch_versions = Channel.empty()

    INPUT_CHECK( params.input )
    .reads
    .set { ch_samples_from_sheet }

    ch_samples_from_sheet
    .branch {
        meta, filename ->
            raw : meta.ext == 'raw'
                return [ meta, filename ]
            mzml :  meta.ext == 'mzml'
                return [ meta, filename ]
            other : true }
        .set { ms_files }

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Input fasta file
    Channel.fromPath( params.fasta )
        .combine( ch_samples_from_sheet )
        .flatMap{ it -> [tuple(it[1],it[0])] }
        .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
        .set { input_fasta }

    // A warning message will be given when the format differs from the '.raw' or '.mzML' extention
    ms_files.other.subscribe { row -> log.warn("Unknown format for entry " + row[3] + " in provided sample sheet, line will be ignored."); exit 1 }

    if ( params.include_proteins_from_vcf ) {
        // Combine the vcf information with the meta information
        ch_vcf = input_fasta
            .map{ it -> [it[0].sample, it[0], it[1]] }
            .join( ch_vcf_from_sheet )
            .map(it -> [it[1], it[2], it[3]])
        // If specified translate variants to proteins and include in reference fasta
        GENERATE_PROTEINS_FROM_VCF( ch_vcf )
        // ch_versions = ch_versions.mix(GENERATE_PROTEINS_FROM_VCF.out.versions.first().ifEmpty(null))
        ch_fasta_file = GENERATE_PROTEINS_FROM_VCF.out.vcf_fasta
    } else {
        ch_fasta_file = input_fasta
    }

    if (!params.skip_decoy_generation) {
        // Generate reversed decoy database
        OPENMS_DECOYDATABASE(ch_fasta_file)
        ch_versions = ch_versions.mix(OPENMS_DECOYDATABASE.out.versions.first().ifEmpty(null))
        ch_decoy_db = OPENMS_DECOYDATABASE.out.decoy
    } else {
        ch_decoy_db = ch_fasta_file
    }

    // Raw file conversion
    OPENMS_THERMORAWFILEPARSER(ms_files.raw)
    ch_versions = ch_versions.mix(OPENMS_THERMORAWFILEPARSER.out.versions.first().ifEmpty(null))

    if ( params.run_centroidisation ) {
        // Optional: Run Peak Picking as Preprocessing
        OPENMS_PEAKPICKERHIRES(ms_files.mzml)
        ch_versions = ch_versions.mix(OPENMS_PEAKPICKERHIRES.out.versions.first().ifEmpty(null))
        ch_mzml_file = OPENMS_PEAKPICKERHIRES.out.mzml
    } else {
        ch_mzml_file = ms_files.mzml
    }

    // Run comet database search
    OPENMS_COMETADAPTER(
        OPENMS_THERMORAWFILEPARSER.out.mzml
                .mix(ch_mzml_file)
                .join(ch_decoy_db, remainder:true))
    ch_versions = ch_versions.mix(OPENMS_COMETADAPTER.out.versions.first().ifEmpty(null))

    // Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(OPENMS_COMETADAPTER.out.idxml.join(ch_decoy_db))
    ch_versions = ch_versions.mix(OPENMS_PEPTIDEINDEXER.out.versions.first().ifEmpty(null))

    if(!params.skip_quantification) {
        // Calculate fdr for id based alignment
        OPENMS_FALSEDISCOVERYRATE(OPENMS_PEPTIDEINDEXER.out.idxml)
        ch_versions = ch_versions.mix(OPENMS_FALSEDISCOVERYRATE.out.versions.first().ifEmpty(null))
        // Filter fdr for id based alignment
        OPENMS_IDFILTER_FOR_ALIGNMENT(OPENMS_FALSEDISCOVERYRATE.out.idxml
            .flatMap { it -> [tuple(it[0], it[1], null)]})
        ch_versions = ch_versions.mix(OPENMS_IDFILTER_FOR_ALIGNMENT.out.versions.first().ifEmpty(null))

        ch_grouped_fdr_filtered = OPENMS_IDFILTER_FOR_ALIGNMENT.out.idxml
            .map {
                meta, raw ->
                    [[id:meta.sample + "_" + meta.condition, sample:meta.sample, condition:meta.condition, ext:meta.ext], raw]
                }
            .groupTuple(by: [0])

        // Compute alignment rt transformatio
        OPENMS_MAPALIGNERIDENTIFICATION(ch_grouped_fdr_filtered)
        ch_versions = ch_versions.mix(OPENMS_MAPALIGNERIDENTIFICATION.out.versions.first().ifEmpty(null))
        // Intermediate step to join RT transformation files with mzml and idxml channels
        ms_files.mzml
        .mix(OPENMS_THERMORAWFILEPARSER.out.mzml)
        .mix(ch_mzml_file)
        .join(
            OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
                .transpose()
                .flatMap {
                    meta, trafoxml ->
                        ident = trafoxml.baseName.split('_-_')[0]
                        [[[id:ident, sample:meta.sample, condition:meta.condition, ext:meta.ext], trafoxml]]
                }, by: [0] )
        .set { joined_trafos_mzmls }

        OPENMS_PEPTIDEINDEXER.out.idxml
        .join(
            OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
                .transpose()
                .flatMap {
                    meta, trafoxml ->
                        ident = trafoxml.baseName.split('_-_')[0]
                        [[[id:ident, sample:meta.sample, condition:meta.condition, ext:meta.ext], trafoxml]]
                }, by: [0] )
        .set { joined_trafos_ids }

        // Align mzML files using trafoXMLs
        OPENMS_MAPRTTRANSFORMERMZML(joined_trafos_mzmls)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERMZML.out.versions.first().ifEmpty(null))
        // Align unfiltered idXMLfiles using trafoXMLs
        OPENMS_MAPRTTRANSFORMERIDXML(joined_trafos_ids)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERIDXML.out.versions.first().ifEmpty(null))
        ch_proceeding_idx = OPENMS_MAPRTTRANSFORMERIDXML.out.aligned
            .map {
                meta, raw ->
                [[id:meta.sample + "_" + meta.condition, sample:meta.sample, condition:meta.condition, ext:meta.ext], raw]
            }
            .groupTuple(by: [0])

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
    ch_versions = ch_versions.mix(OPENMS_IDMERGER.out.versions.first().ifEmpty(null))
    // Extract PSM features for Percolator
    OPENMS_PSMFEATUREEXTRACTOR(OPENMS_IDMERGER.out.idxml)
    ch_versions = ch_versions.mix(OPENMS_PSMFEATUREEXTRACTOR.out.versions.first().ifEmpty(null))
    // Run Percolator
    if (params.description_correct_features > 0 && params.klammer) {
        OPENMS_PERCOLATORADAPTER_KLAMMER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER_KLAMMER.out.versions.first().ifEmpty(null))
        ch_percolator_adapter_outcome = OPENMS_PERCOLATORADAPTER_KLAMMER.out.idxml

    } else {
        OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions.first().ifEmpty(null))
        ch_percolator_adapter_outcome = OPENMS_PERCOLATORADAPTER.out.idxml
    }

    // Filter by percolator q-value
    OPENMS_IDFILTER_Q_VALUE(ch_percolator_adapter_outcome.flatMap { it -> [tuple(it[0], it[1], null)]})
    ch_versions = ch_versions.mix(OPENMS_IDFILTER_Q_VALUE.out.versions.first().ifEmpty(null))

    // Refine_fdr_on_predicted_subset
    if ( params.refine_fdr_on_predicted_subset && params.predict_class_1 ) {
        // Run the following subworkflow
        REFINE_FDR_ON_PREDICTED_SUBSET (
            OPENMS_IDFILTER_Q_VALUE.out.idxml,
            OPENMS_PSMFEATUREEXTRACTOR.out.idxml,
            peptides_class_1_alleles
        )
        ch_versions = ch_versions.mix(REFINE_FDR_ON_PREDICTED_SUBSET.out.versions.first().ifEmpty(null))

        // Define the outcome of the paramer to a fixed variable
        filter_q_value = REFINE_FDR_ON_PREDICTED_SUBSET.out.filter_refined_q_value.flatMap { it -> [ tuple(it[0].sample, it[0], it[1]) ] }
    } else {
        // Make sure that the columns that consists of the ID's, sample names and the idXML file names are returned
        filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.idxml.map{ it -> [it[0].sample, it[0], it[1]] }
    }

    if ( !params.skip_quantification) {
        // Combining the necessary information into one channel
        OPENMS_IDFILTER_FOR_ALIGNMENT.out[0]
            .join( OPENMS_MAPRTTRANSFORMERMZML.out[0], by: [0] )
            .map { it -> [it[0].sample, it[0], it[1], it[2]] }
            .combine( filter_q_value , by: [0] )
            .map { it -> [it[1], it[2], it[3], it[5]] }
            .set{ joined_mzmls_ids_quant }
        // Quantify identifications using targeted feature extraction
        OPENMS_FEATUREFINDERIDENTIFICATION(joined_mzmls_ids_quant)
        ch_versions = ch_versions.mix(OPENMS_FEATUREFINDERIDENTIFICATION.out.versions.first().ifEmpty(null))
        // Link extracted features
        OPENMS_FEATURELINKERUNLABELEDKD(
            OPENMS_FEATUREFINDERIDENTIFICATION.out.featurexml
                .flatMap {
                    meta, raw ->
                        [[[id:meta.sample + "_" + meta.condition, sample:meta.sample, condition:meta.condition, ext:meta.ext], raw]]
                }
                .groupTuple(by:[0]))
        ch_versions = ch_versions.mix(OPENMS_FEATURELINKERUNLABELEDKD.out.versions.first().ifEmpty(null))
        // Resolve conflicting ids matching to the same feature
        OPENMS_IDCONFLICTRESOLVER(OPENMS_FEATURELINKERUNLABELEDKD.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_IDCONFLICTRESOLVER.out.versions.first().ifEmpty(null))
        // Assign the outcome of the id conflict resolver as export content
        export_content = OPENMS_IDCONFLICTRESOLVER.out.consensusxml
    } else {
        // Assign the outcome of the filter q value as export content
        export_content = filter_q_value.map { it -> [it[1], it[2]] }
    }

    // Export all information as text to csv
    OPENMS_TEXTEXPORTER(export_content)
    // Export all information as mzTab
    OPENMS_MZTABEXPORTER(export_content)
    //////////////////////////////////////////////////////////////////////////////////////////////
    //  TODO: Replacement of custom scripts with epytope
    ch_predicted_possible_neoepitopes = Channel.empty()
    if ( params.predict_class_1  & !params.skip_quantification ) {
        // If specified predict peptides using MHCFlurry
        PREDICT_PEPTIDES_MHCFLURRY_CLASS_1(
            OPENMS_MZTABEXPORTER.out.mztab
                .map{ it -> [it[0].sample, it[0], it[1]] }
                .combine( peptides_class_1_alleles, by:0)
                .map( it -> [it[1], it[2], it[3]])
            )
        ch_versions = ch_versions.mix(PREDICT_PEPTIDES_MHCFLURRY_CLASS_1.out.versions.first().ifEmpty(null))
        if ( params.include_proteins_from_vcf ) {
            // Predict all possible neoepitopes from vcf
            PREDICT_POSSIBLE_NEOEPITOPES(peptides_class_1_alleles.join(ch_vcf_from_sheet, by:0, remainder:true))
            ch_versions = ch_versions.mix(PREDICT_POSSIBLE_NEOEPITOPES.out.versions.first().ifEmpty(null))
            ch_predicted_possible_neoepitopes = PREDICT_POSSIBLE_NEOEPITOPES.out.csv
            // Resolve found neoepitopes
            RESOLVE_FOUND_NEOEPITOPES(
                OPENMS_MZTABEXPORTER.out.mztab
                    .map{ it -> [it[0].sample, it[0], it[1]] }
                    .combine( ch_predicted_possible_neoepitopes, by:0, remainder:true)
                    .map( it -> [it[1], it[2], it[3]])
                )
            ch_versions = ch_versions.mix(RESOLVE_FOUND_NEOEPITOPES.out.versions.first().ifEmpty(null))
            // Predict class 1 neoepitopes MHCFlurry
            PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1(peptides_class_1_alleles.join(RESOLVE_FOUND_NEOEPITOPES.out.csv, by:0))
            ch_versions = ch_versions.mix(PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1.out.versions.first().ifEmpty(null))
        }
    }

    ch_predicted_possible_neoepitopes_II = Channel.empty()
    if ( params.predict_class_2  & !params.skip_quantification ) {
        // Preprocess found peptides for MHCNuggets prediction class 2
        PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2(OPENMS_MZTABEXPORTER.out.mztab)
        ch_versions = ch_versions.mix(PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out.versions.first().ifEmpty(null))
        // Predict found peptides using MHCNuggets class 2
        PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2(
            PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out.preprocessed
                .map{ it -> [it[0].sample, it[0], it[1]] }
                .join(peptides_class_2_alleles, by:0)
                .map( it -> [it[1], it[2], it[3]])
        )
        ch_versions = ch_versions.mix(PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2.out.versions.first().ifEmpty(null))
        // Postprocess predicted MHCNuggets peptides class 2
        POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2( PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2.out.csv.join(PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out.geneID, by:0) )
        ch_versions = ch_versions.mix(POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out.versions.first().ifEmpty(null))

        if ( params.include_proteins_from_vcf ) {
            // Predict all possible class 2 neoepitopes from vcf
            PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES(peptides_class_2_alleles.join(ch_vcf_from_sheet, by:0, remainder:true))
            ch_versions = ch_versions.mix(PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES.out.versions.first().ifEmpty(null))
            ch_predicted_possible_neoepitopes_II = PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES.out.csv
            // Resolve found class 2 neoepitopes
            RESOLVE_FOUND_CLASS_2_NEOEPITOPES(
                OPENMS_MZTABEXPORTER.out.mztab
                    .map{ it -> [it[0].sample, it[1]] }
                    .combine( ch_predicted_possible_neoepitopes_II, by:0, remainder:true)
            )
            ch_versions = ch_versions.mix(RESOLVE_FOUND_CLASS_2_NEOEPITOPES.out.versions.first().ifEmpty(null))
            // Preprocess resolved neoepitopes in a format that MHCNuggets understands
            PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2(RESOLVE_FOUND_CLASS_2_NEOEPITOPES.out.csv)
            ch_versions = ch_versions.mix(PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.versions.first().ifEmpty(null))
            // Predict class 2 MHCNuggets
            PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2(PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.preprocessed.join(peptides_class_2_alleles, by:0))
            ch_versions = ch_versions.mix(PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.versions.first().ifEmpty(null))
            // Class 2 MHCNuggets Postprocessing
            POSTPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2(RESOLVE_FOUND_CLASS_2_NEOEPITOPES.out.csv.join(PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.csv, by:0))
            ch_versions = ch_versions.mix(POSTPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.versions.first().ifEmpty(null))
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////
    if ( params.predict_RT ) {
        filter_q_value = filter_q_value.map{ it -> [it[1], it[2]] }
        // Train Retention Times Predictor
        OPENMS_RTMODEL(filter_q_value)
        ch_versions = ch_versions.mix(OPENMS_RTMODEL.out.versions.first().ifEmpty(null))
        // Retention Times Predictor Found Peptides
        OPENMS_RTPREDICT_FOUND_PEPTIDES(filter_q_value.join(OPENMS_RTMODEL.out.complete, by:[0]))
        ch_versions = ch_versions.mix(OPENMS_RTPREDICT_FOUND_PEPTIDES.out.versions.first().ifEmpty(null))
        // Retention Times Predictor possible Neoepitopes
        OPENMS_RTPREDICT_NEOEPITOPES(ch_predicted_possible_neoepitopes.mix(ch_predicted_possible_neoepitopes_II).join(OPENMS_RTMODEL.out.complete, by:[0]))
        ch_versions = ch_versions.mix(OPENMS_RTPREDICT_FOUND_PEPTIDES.out.versions.first().ifEmpty(null))
    }

    //
    // MODULE: Pipeline reporting
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile()
    )

}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////
workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    NfcoreTemplate.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
