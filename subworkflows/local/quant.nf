/*
 * Perform the quantification of the samples when the parameter --skip_quantification is not provided
 * This workflow splits the merged percolator output into the individual runs and filters them based on the q-value
 * It then aligns the retention times of the runs and merges the idxml files together to use them as id_ext in featurefinder
 * Finally, it performs the quantification and emits the consensusXML file
 */
include { OPENMS_IDRIPPER }                                                 from '../../modules/local/openms_idripper'
include { OPENMS_IDSCORESWITCHER }                                          from '../../modules/local/openms_idscoreswitcher'
include { PYOPENMS_IDFILTER }                                               from '../../modules/local/pyopenms_idfilter'
include { OPENMS_IDMERGER as OPENMS_IDMERGER_QUANT }                        from '../../modules/local/openms_idmerger'

include { MAP_ALIGNMENT }                                                   from './map_alignment'
include { PROCESS_FEATURE }                                                 from './process_feature'

// Sort closure for merging and splitting files
def sortById = { a, b -> a.id <=> b.id }

workflow QUANT {
    take:
        merge_meta_map
        merged_pout
        filter_q_value
        mzml

    main:
        ch_versions = Channel.empty()
        // Rip post-percolator idXML files and manipulate such that we end up with [meta_run1, idxml_run1, pout_filtered] [meta_run2, idxml_run2, pout_filtered] ...
        OPENMS_IDRIPPER( merged_pout ).ripped
                .join( merge_meta_map )
                .join( filter_q_value )
                .map { group_meta, ripped, meta, fdrfiltered -> [meta, ripped, fdrfiltered] }
                .transpose()
                .set { ch_ripped_pout }
        ch_versions = ch_versions.mix(OPENMS_IDRIPPER.out.versions.ifEmpty(null))

        // Switch to xcorr for filtering since q-values are set to 1 with peptide-level-fdr
        if (params.fdr_level == 'peptide_level_fdrs'){
            ch_runs_to_be_filtered = OPENMS_IDSCORESWITCHER( ch_ripped_pout ).switched_idxml
            ch_versions = ch_versions.mix(OPENMS_IDSCORESWITCHER.out.versions.ifEmpty(null))
        } else {
            ch_runs_to_be_filtered = ch_ripped_pout
        }

        // Filter runs based on fdr filtered coprocessed percolator output.
        // NOTE: This is an alternative filtering method that will be replaced by IDFilter with new release of OpenMS
        PYOPENMS_IDFILTER( ch_runs_to_be_filtered ).filtered
                .map { meta, idxml -> [[id:meta.sample + '_' + meta.condition], [id:meta.id, file:idxml]] }
                .groupTuple( sort: sortById )
                .map { meta, idxml -> [meta, idxml.file] }
                .set { ch_runs_to_be_aligned }
        ch_versions = ch_versions.mix(PYOPENMS_IDFILTER.out.versions.ifEmpty(null))

        // Align retention times of runs
        MAP_ALIGNMENT(
            ch_runs_to_be_aligned,
            mzml,
            merge_meta_map
        )
        ch_versions = ch_versions.mix( MAP_ALIGNMENT.out.versions.ifEmpty(null) )

        // We need to merge groupwise the aligned idxml files together to use them as id_ext in featurefinder
        OPENMS_IDMERGER_QUANT( MAP_ALIGNMENT.out.aligned_idxml
                                    .map { meta, aligned_idxml -> [[id: meta.sample + '_' + meta.condition], aligned_idxml] }
                                    .groupTuple())
        ch_versions = ch_versions.mix(OPENMS_IDMERGER_QUANT.out.versions.ifEmpty(null))

        // Manipulate channels such that we end up with : [meta, mzml, run_idxml, merged_runs_idxml]
        MAP_ALIGNMENT.out.aligned_mzml
                .join( MAP_ALIGNMENT.out.aligned_idxml )
                .map { meta, mzml, idxml -> [[id: meta.sample + '_' + meta.condition], meta, [id:meta.id, file:mzml], [id:meta.id, file:idxml]] }
                .groupTuple( sort: sortById )
                .map { group_meta, meta, mzml, idxml -> [group_meta, meta, mzml.file, idxml.file] }
                .join( OPENMS_IDMERGER_QUANT.out.idxml )
                .map { group_meta, meta, mzml, idxml, merged_idxml -> [meta, mzml, idxml, merged_idxml] }
                .transpose()
                .set { ch_runs_to_be_quantified }

        PROCESS_FEATURE ( ch_runs_to_be_quantified )
        ch_versions = ch_versions.mix(PROCESS_FEATURE.out.versions.ifEmpty(null))

    emit:
        consensusxml = PROCESS_FEATURE.out.consensusxml
        versions = ch_versions
}
