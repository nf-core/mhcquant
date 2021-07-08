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

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////
// Input/outpout options
if (params.input)   { sample_sheet = file(params.input) }   else { exit 1, 'Input samplesheet was not specified!' }
if (params.fasta)   { params.fasta = params.fasta }         else { exit 1, 'No database fasta was provided, make sure you have used the '--fasta' option.' }
if (params.outdir)  { params.outdir  = './results' }        else { log.warn 'Results into \'./results\'.\nIf you want to define a result directory, please use the --outdir option.' }


// Read the content of the sample sheet
// Reading to content of the sample sheet and defining the different columns for further use
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
                .flatMap {it -> [tuple("id", it[0].toString(), it[1])] }
                .set { peptides_class_1_alleles }
        }

        // Allele class 2
        if( params.predict_class_2){
            ch_alleles_from_sheet.classII
                .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
                .flatMap { it -> [tuple("id", it[0].toString(), it[1])] }
                .set { peptides_class_2_alleles }
        }
}

// Variant 
if (params.include_proteins_from_vcf)  {
    // TODO: test with a VCF sheet
    vcf_sheet = file(params.vcf_sheet, checkIfExists: true)
    Channel.from( vcf_sheet )
        .splitCsv(header: ['Sample', 'VCF_FileName'], sep:'\t', skip: 1)
        .map { col -> tuple("${col.Sample}", file("${col.VCF_FileName}"),) }
        .set { ch_vcf_from_sheet }

    ch_vcf_from_sheet
        .ifEmpty { exit 1, "params.vcf_sheet was empty - no vcf input file supplied" }
        .flatMap {it -> [tuple("id", it[0].toString(), it[1])] }
        .set { input_vcf }
}

// Database options
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

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////
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
id_filter_qvalue_options.suffix = "_all_ids_merged_psm_perc_filtered"





include { hasExtension }                                    from '../modules/local/functions'

// include { INPUT_CHECK }                                     from '../modules/local/subworkflow/input_check'                                  addParams( options: [:] )           
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

include { REFINE_FDR_ON_PREDICTED_SUBSET }                  from '../modules/local/subworkflow/refine_fdr_on_predicted_subset'               addParams( run_percolator_options : percolator_adapter_options, filter_options: id_filter_options, whitelist_filter_options: id_filter_whitelist_options)

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

include { GET_SOFTWARE_VERSIONS }                           from '../modules/local/get_software_versions'                            addParams( options: [publish_files : ['csv':'']]   )

////////////////////////////////////////////////////
/* --              CREATE CHANNELS             -- */
////////////////////////////////////////////////////
// Input fasta file
Channel.fromPath( params.fasta )
    .combine( ch_samples_from_sheet )
    .flatMap{it -> [tuple(it[1].toInteger(),it[2],it[0])]}      //maps tuple to val("id"), val("Sample"), val("Condition"), val("ReplicateFileName")
    .ifEmpty { exit 1, "params.fasta was empty - no input file supplied" }
    .set { input_fasta }

params.summary_params = [:]

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow MHCQUANT {

    //INPUT_CHECK(sample_sheet)
    //.map { meta, raw -> meta.id = meta.id.split('_')[0..-1].join('_')
    //[ meta, raw ] }
    //.groupTuple(by: [0])
    //.map { it -> [ it[0], it[1].flatten() ] }
    //.set { ch_samples_from_sheet }

    ch_software_versions = Channel.empty()
    // A warning message will be given when the format differs from the '.raw' or '.mzML' extention
    ms_files.other.subscribe { row -> log.warn("Unknown format for entry " + row[3] + " in provided sample sheet, line will be ignored."); exit 1 }

    // ////////////////////////////////////////////////////////////////////////////////////////////
    if ( params.include_proteins_from_vcf ) {
        // If specified translate variants to proteins and include in reference fasta
        GENERATE_PROTEINS_FROM_VCF(input_fasta.combine(input_vcf, by:1))
        ch_fasta_file = GENERATE_PROTEINS_FROM_VCF.out.vcf_fasta
    } else {
        ch_fasta_file = input_fasta
    }
    
    if (!params.skip_decoy_generation) {
        // Generate reversed decoy database
        OPENMS_DECOYDATABASE(ch_fasta_file)
        ch_decoy_db = OPENMS_DECOYDATABASE.out.decoy
    } else {
        ch_decoy_db = ch_fasta_file.flatMap { it -> [tuple(it[0], it[1], null)]}
    }

    // Raw file conversion
    OPENMS_THERMORAWFILEPARSER(ms_files.raw)
    ch_software_versions = ch_software_versions.mix(OPENMS_THERMORAWFILEPARSER.out.version.first().ifEmpty(null))

    if ( params.run_centroidisation ) {
        // Optional: Run Peak Picking as Preprocessing
        OPENMS_PEAKPICKERHIRES(ms_files.mzML)
        ch_mzml_file = OPENMS_PEAKPICKERHIRES.out.mzml
    } else {
        ch_mzml_file = ms_files.mzML
    }
    
    // Run comet database search
    OPENMS_COMETADAPTER(
        OPENMS_THERMORAWFILEPARSER.out.mzml
                .mix(ch_mzml_file)
                .join(ch_decoy_db, by:1, remainder:true))
    ch_software_versions = ch_software_versions.mix(OPENMS_COMETADAPTER.out.version.first().ifEmpty(null))

    // Index decoy and target hits
    OPENMS_PEPTIDEINDEXER(OPENMS_COMETADAPTER.out.idxml.join(ch_decoy_db, by:1))

    if(!params.skip_quantification) {
        // Calculate fdr for id based alignment
        OPENMS_FALSEDISCOVERYRATE(OPENMS_PEPTIDEINDEXER.out[0])
        // Filter fdr for id based alignment
        OPENMS_IDFILTER_FOR_ALIGNMENT(OPENMS_FALSEDISCOVERYRATE.out.idxml
        .flatMap { it -> [tuple(it[0], it[1], it[2], it[3], null)]}) 
        // Compute alignment rt transformation
        OPENMS_MAPALIGNERIDENTIFICATION(OPENMS_IDFILTER_FOR_ALIGNMENT.out.idxml
                .groupTuple(by: 1))

        // Intermediate step to join RT transformation files with mzml and idxml channels
        ms_files.mzML
        .mix(OPENMS_THERMORAWFILEPARSER.out.mzml)
        .mix(ch_mzml_file)
        .flatMap { it -> [tuple(it[0].toInteger(), it[1], it[2], it[3])] }
        .join(OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml.transpose().flatMap{ it -> [tuple(it[1].baseName.split('_-_')[0].toInteger(), it[0], it[1])] }, by: [0,1])
        .set { joined_trafos_mzmls }

        OPENMS_PEPTIDEINDEXER.out.idxml
        .flatMap { it -> [tuple(it[0].toInteger(), it[1], it[2], it[3])] }
        .join(OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml.transpose().flatMap{ it -> [tuple(it[1].baseName.split('_-_')[0].toInteger(), it[0], it[1])] }, by: [0,1])
        .set { joined_trafos_ids }

        // Align mzML files using trafoXMLs
        OPENMS_MAPRTTRANSFORMERMZML(joined_trafos_mzmls)
        // Align unfiltered idXMLfiles using trafoXMLs
        OPENMS_MAPRTTRANSFORMERIDXML(joined_trafos_ids)

        ch_proceeding_idx = OPENMS_MAPRTTRANSFORMERIDXML.out[0].groupTuple(by: 1)
    } else {
        ch_proceeding_idx = OPENMS_PEPTIDEINDEXER.out.idxml.groupTuple(by: 1)
    }

    // Merge aligned idXMLfiles
    OPENMS_IDMERGER(ch_proceeding_idx)
    // Extract PSM features for Percolator
    OPENMS_PSMFEATUREEXTRACTOR(OPENMS_IDMERGER.out.idxml)
    // Run Percolator
    if (params.description_correct_features > 0 && params.klammer) {
        OPENMS_PERCOLATORADAPTER_KLAMMER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_percolator_adapter_outcome = OPENMS_PERCOLATORADAPTER_KLAMMER.out[0]

    } else {
        OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_percolator_adapter_outcome = OPENMS_PERCOLATORADAPTER.out[0]
    }
    
    // Filter by percolator q-value
    OPENMS_IDFILTER_Q_VALUE(ch_percolator_adapter_outcome.flatMap { it -> [tuple(it[0], it[1], it[2], it[3], null)]}) 
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Refine_fdr_on_predicted_subset
    if ( params.refine_fdr_on_predicted_subset) {
        // Run the following subworkflow
        REFINE_FDR_ON_PREDICTED_SUBSET (
            OPENMS_IDFILTER_Q_VALUE.out[0],
            OPENMS_PSMFEATUREEXTRACTOR.out[0],
            peptides_class_1_alleles
        )

        // Define the outcome of the paramer to a fixed variable
        filter_q_value = REFINE_FDR_ON_PREDICTED_SUBSET.out[0].flatMap { it -> [ tuple( it[0], it[1], it[3] ) ] }
    } else {
        // Make sure that the columns that consists of the ID's, sample names and the idXML file names are returned
        filter_q_value = OPENMS_IDFILTER_Q_VALUE.out[0]
        .flatMap { it -> [ tuple( it[0], it[1], it[3] ) ] }
    }

    OPENMS_IDFILTER_FOR_ALIGNMENT.out[0]
    .flatMap { it -> [ tuple( it[0], it[1], it[2], it[3] ) ] }
    .join( OPENMS_MAPRTTRANSFORMERMZML.out[0], by: [0,1,2] )
    .combine( filter_q_value , by: 1 )
    .set{ joined_mzmls_ids_quant }

    if ( !params.skip_quantification) {
        // Quantify identifications using targeted feature extraction
        OPENMS_FEATUREFINDERIDENTIFICATION(joined_mzmls_ids_quant)
        // Link extracted features
        OPENMS_FEATURELINKERUNLABELEDKD(OPENMS_FEATUREFINDERIDENTIFICATION.out.featurexml.groupTuple(by:1))  
        // Resolve conflicting ids matching to the same feature
        OPENMS_IDCONFLICTRESOLVER(OPENMS_FEATURELINKERUNLABELEDKD.out.consensusxml)
    }

    // Export all information as text to csv
    OPENMS_TEXTEXPORTER(OPENMS_IDCONFLICTRESOLVER.out.consensusxml)
    // Export all information as mzTab
    OPENMS_MZTABEXPORTER(
        OPENMS_IDCONFLICTRESOLVER.out.consensusxml
        .flatMap { it -> [tuple("id", it[0], it[0], it[1])]}
    )

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  TODO: Replacement of custom scripts with epytope
    ch_predicted_possible_neoepitopes = Channel.empty()
    if ( params.predict_class_1 ) {
        // If specified predict peptides using MHCFlurry
        PREDICT_PEPTIDES_MHCFLURRY_CLASS_1(OPENMS_MZTABEXPORTER.out.mztab.combine(peptides_class_1_alleles, by:1))
        ch_software_versions = ch_software_versions.mix(PREDICT_PEPTIDES_MHCFLURRY_CLASS_1.out.version.first().ifEmpty(null))
        if ( params.include_proteins_from_vcf ) {
            // Predict class 1 neoepitopes MHCFlurry
            PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1(peptides_class_1_alleles.join(RESOLVE_FOUND_NEOEPITOPES.out.csv, by:1))
            // Predict all possible neoepitopes from vcf
            PREDICT_POSSIBLE_NEOEPITOPES(peptides_class_1_alleles.join(input_vcf, by:[0,1], remainder:true))
            ch_predicted_possible_neoepitopes = PREDICT_POSSIBLE_NEOEPITOPES.out[0]
            // Resolve found neoepitopes
            RESOLVE_FOUND_NEOEPITOPES(OPENMS_MZTABEXPORTER.out.mztab.join(ch_predicted_possible_neoepitopes, by:[0,1], remainder:true))
        }
    }

    ch_predicted_possible_neoepitopes_II = Channel.empty()
    if ( params.predict_class_2 ) {
        // Preprocess found peptides for MHCNuggets prediction class 2
        PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2(OPENMS_MZTABEXPORTER.out.mztab)
        // Predict found peptides using MHCNuggets class 2
        PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2(PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out[0].join(peptides_class_2_alleles, by:1))
        // Postprocess predicted MHCNuggets peptides class 2
        POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2(PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2.out.csv.join(PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2.out[1], by:1))

        if ( params.include_proteins_from_vcf ) { 
            // Predict all possible class 2 neoepitopes from vcf
            PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES(peptides_class_2_alleles.join(input_vcf, by:[0,1], remainder:true))            
            ch_predicted_possible_neoepitopes_II = PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES.out[0]
            // Resolve found class 2 neoepitopes
            RESOLVE_FOUND_CLASS_2_NEOEPITOPES(OPENMS_MZTABEXPORTER.out.mztab.join(ch_predicted_possible_neoepitopes_II, by:[0,1], remainder:true)) 
            // Preprocess resolved neoepitopes in a format that MHCNuggets understands
            PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2(RESOLVE_FOUND_CLASS_2_NEOEPITOPES.out.csv)
            // Predict class 2 MHCNuggets
            PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2(PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.preprocessed.join(peptides_class_2_alleles, by:1))
            // Class 2 MHCNuggets Postprocessing
            POSTPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2(RESOLVE_FOUND_CLASS_2_NEOEPITOPES.out.csv.join(PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2.out.csv, by:1))
            // If there was no prediction performed on class 1 
            if ( !params.predict_class_1 ) {
                // Add information to software versions
                ch_software_versions = ch_software_versions.mix(PREDICT_POSSIBLE_CLASS_2_NEOEPITOPES.out.version.ifEmpty(null))
            }
        } else if (!params.include_proteins_from_vcf && !params.predict_class_1) {
            // Add the information to software versions
            ch_software_versions = ch_software_versions.mix(PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2.out.version.first().ifEmpty(null))
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////

    if ( params.predict_RT ) {        
        // Train Retention Times Predictor
        OPENMS_RTMODEL(filter_q_value)
        // Retention Times Predictor Found Peptides
        OPENMS_RTPREDICT_FOUND_PEPTIDES(filter_q_value.join(OPENMS_RTMODEL.out.complete, by:[0,1]))  
        // Retention Times Predictor possible Neoepitopes
        OPENMS_RTPREDICT_NEOEPITOPES(ch_predicted_possible_neoepitopes.mix(ch_predicted_possible_neoepitopes_II).join(OPENMS_RTMODEL.out.complete, by:[0,1])) 
    }

    /*
     * MODULE: Pipeline reporting
     */
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

}


////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////
workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
