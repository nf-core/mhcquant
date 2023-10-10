/*
 * Perform the class 2 prediction when the parameter --predict_class_2 is provided and --skip_quantification is not
 */

include { MHCNUGGETS_PEPTIDESCLASS2PRE }                                        from '../../modules/local/mhcnuggets_peptidesclass2pre'
include { MHCNUGGETS_PREDICTPEPTIDESCLASS2 }                                    from '../../modules/local/mhcnuggets_predictpeptidesclass2'
include { MHCNUGGETS_PEPTIDESCLASS2POST }                                       from '../../modules/local/mhcnuggets_peptidesclass2post'
include { PREDICT_POSSIBLE_NEOEPITOPES as PREDICT_POSSIBLE_CLASS2_NEOEPITOPES } from '../../modules/local/predict_possible_neoepitopes'
include { RESOLVE_FOUND_NEOEPITOPES as RESOLVE_FOUND_CLASS2_NEOEPITOPES }       from '../../modules/local/resolve_found_neoepitopes'
include { MHCNUGGETS_NEOEPITOPESCLASS2PRE }                                     from '../../modules/local/mhcnuggets_neoepitopesclass2pre'
include { MHCNUGGETS_PREDICTNEOEPITOPESCLASS2 }                                 from '../../modules/local/mhcnuggets_predictneoepitopesclass2'
include { MHCNUGGETS_NEOEPITOPESCLASS2POST }                                    from '../../modules/local/mhcnuggets_neoepitopesclass2post'

workflow PREDICT_CLASS2 {
    take:
        mztab
        peptides_class_2_alleles
        ch_vcf_from_sheet

    main:
        ch_versions = Channel.empty()
        ch_predicted_possible_neoepitopes = Channel.empty()
        alleles = peptides_class_2_alleles.map{meta, alleles -> [[id:meta], alleles]}

        // Preprocess found peptides for MHCNuggets prediction class 2
        MHCNUGGETS_PEPTIDESCLASS2PRE(mztab)
        ch_versions = ch_versions.mix(MHCNUGGETS_PEPTIDESCLASS2PRE.out.versions.first().ifEmpty(null))

        // Predict found peptides using MHCNuggets class 2
        MHCNUGGETS_PREDICTPEPTIDESCLASS2(
            MHCNUGGETS_PEPTIDESCLASS2PRE.out.preprocessed
                .join(alleles)
        )
        ch_versions = ch_versions.mix(MHCNUGGETS_PREDICTPEPTIDESCLASS2.out.versions.first().ifEmpty(null))
        // Postprocess predicted MHCNuggets peptides class 2
        MHCNUGGETS_PEPTIDESCLASS2POST( MHCNUGGETS_PREDICTPEPTIDESCLASS2.out.csv.join(MHCNUGGETS_PEPTIDESCLASS2PRE.out.geneID, by:0) )
        ch_versions = ch_versions.mix(MHCNUGGETS_PEPTIDESCLASS2POST.out.versions.first().ifEmpty(null))
        if ( params.include_proteins_from_vcf ) {
            // Predict all possible class 2 neoepitopes from vcf
            PREDICT_POSSIBLE_CLASS2_NEOEPITOPES(alleles.combine(ch_vcf_from_sheet, by:0))
            ch_versions = ch_versions.mix(PREDICT_POSSIBLE_CLASS2_NEOEPITOPES.out.versions.first().ifEmpty(null))
            ch_predicted_possible_neoepitopes = PREDICT_POSSIBLE_CLASS2_NEOEPITOPES.out.csv
            // Resolve found class 2 neoepitopes
            RESOLVE_FOUND_CLASS2_NEOEPITOPES(
                mztab
                    .map{ it -> [it[0].sample, it[1]] }
                    .combine( ch_predicted_possible_neoepitopes, by:0)
            )
            ch_versions = ch_versions.mix(RESOLVE_FOUND_CLASS2_NEOEPITOPES.out.versions.first().ifEmpty(null))
            // Preprocess resolved neoepitopes in a format that MHCNuggets understands
            MHCNUGGETS_NEOEPITOPESCLASS2PRE(RESOLVE_FOUND_CLASS2_NEOEPITOPES.out.csv)
            ch_versions = ch_versions.mix(MHCNUGGETS_NEOEPITOPESCLASS2PRE.out.versions.first().ifEmpty(null))
            // Predict class 2 MHCNuggets
            MHCNUGGETS_PREDICTNEOEPITOPESCLASS2(MHCNUGGETS_NEOEPITOPESCLASS2PRE.out.preprocessed.join(alleles, by:0))
            ch_versions = ch_versions.mix(MHCNUGGETS_PREDICTNEOEPITOPESCLASS2.out.versions.first().ifEmpty(null))
            // Class 2 MHCNuggets Postprocessing
            MHCNUGGETS_NEOEPITOPESCLASS2POST(RESOLVE_FOUND_CLASS2_NEOEPITOPES.out.csv.join(MHCNUGGETS_PREDICTNEOEPITOPESCLASS2.out.csv, by:0))
            ch_versions = ch_versions.mix(MHCNUGGETS_NEOEPITOPESCLASS2POST.out.versions.first().ifEmpty(null))
        }

    emit:
        // Define the information that is returned by this workflow
        versions = ch_versions
        ch_predicted_possible_neoepitopes = ch_predicted_possible_neoepitopes
}
