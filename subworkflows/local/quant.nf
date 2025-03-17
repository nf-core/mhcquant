/*
 * Perform the quantification of the samples when the parameter --quantify is provided
 * This workflow splits the merged percolator output into the individual runs and filters them based on the q-value
 * It then aligns the retention times of the runs and merges the idxml files together to use them as id_ext in featurefinder
 * Finally, it performs the quantification and emits the consensusXML file
 */

include { OPENMS_IDRIPPER                          } from '../../modules/nf-core/openms/idripper/main'
include { OPENMS_IDSCORESWITCHER                   } from '../../modules/nf-core/openms/idscoreswitcher/main'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_QUANT } from '../../modules/nf-core/openms/idfilter/main'
include { OPENMS_IDMERGER as OPENMS_IDMERGER_QUANT } from '../../modules/nf-core/openms/idmerger/main'
include { OPENMS_MZTABEXPORTER                     } from '../../modules/local/openms_mztabexporter'

include { MAP_ALIGNMENT                            } from './map_alignment'
include { PROCESS_FEATURE                          } from './process_feature'

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
        // Split post-percolator idXML files and manipulate such that we end up with [meta_run1, idxml_run1] [meta_run2, idxml_run2] ...
        // We need to make sure that the order of the runs is the same as in the mzml files since IDRipper always sorts the runs
        // (and nextflow does not guarantee the order of the maps in merged_meta_map)
        OPENMS_IDRIPPER( merged_pout ).idxmls
                .flatMap { group_meta, idxmls -> idxmls.collect { idxml -> [[spectra: idxml.baseName], idxml] } }
                // join on file basename to make sure that the order of the runs is the same as in the mzml files
                // Is there a smoother way to do this?
                .join( merge_meta_map
                        .flatMap { group_meta, metas -> metas }
                        .map { meta -> [[spectra:meta.spectra], meta]} )
                .map { spectra, idxmls, meta -> [meta, idxmls] }
                .set { ch_ripped_idxml }
        ch_versions = ch_versions.mix(OPENMS_IDRIPPER.out.versions)

        // Switch to xcorr for filtering since q-values are set to 1 with peptide-level-fdr
        if (params.fdr_level == 'peptide_level_fdrs'){
            ch_runs_score_switched = OPENMS_IDSCORESWITCHER( ch_ripped_idxml ).idxml
            ch_versions = ch_versions.mix(OPENMS_IDSCORESWITCHER.out.versions)
        } else {
            ch_runs_score_switched = ch_ripped_idxml
        }

        // Manipulate such that [meta_run1, idxml_run1, pout_group1], [meta_run2, idxml_run2, pout_group1] ...
        ch_runs_score_switched
            // Nextflow can only combine/join on the exact groupKey object, merge_id is not sufficient
            .map { meta, idxml -> [groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count) , meta, idxml] }
            .combine(filter_q_value, by:0)
            .map { group_meta, meta, idxml, q_value -> [meta, idxml, q_value] }
            .set { ch_runs_to_filter}

        // Filter runs based on fdr filtered coprocessed percolator output.
        OPENMS_IDFILTER_QUANT( ch_runs_to_filter ).filtered
                .map { meta, idxml -> [ groupKey([id:"${meta.sample}_${meta.condition}"], meta.group_count), idxml] }
                .groupTuple()
                .set { ch_runs_to_be_aligned }
        ch_versions = ch_versions.mix(OPENMS_IDFILTER_QUANT.out.versions)

        // Align retention times of runs
        MAP_ALIGNMENT(
            ch_runs_to_be_aligned,
            mzml,
            merge_meta_map
        )
        ch_versions = ch_versions.mix( MAP_ALIGNMENT.out.versions )

        // We need to merge groupwise the aligned idxml files together to use them as id_ext in featurefinder
        OPENMS_IDMERGER_QUANT( MAP_ALIGNMENT.out.aligned_idxml
                                    .map { meta, aligned_idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), aligned_idxml] }
                                    .groupTuple())
        ch_versions = ch_versions.mix(OPENMS_IDMERGER_QUANT.out.versions)

        // Manipulate channels such that we end up with : [meta, mzml, run_idxml, merged_runs_idxml]
        MAP_ALIGNMENT.out.aligned_mzml
                .join( MAP_ALIGNMENT.out.aligned_idxml )
                .map { meta, mzml, idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), meta, mzml, idxml] }
                .groupTuple()
                .map { group_meta, meta, mzml, idxml -> [group_meta, meta, mzml, idxml] }
                .join( OPENMS_IDMERGER_QUANT.out.idxml )
                .map { group_meta, meta, mzml, idxml, merged_idxml -> [meta, mzml, idxml, merged_idxml] }
                .transpose()
                .set { ch_runs_to_be_quantified }

        PROCESS_FEATURE ( ch_runs_to_be_quantified )
        ch_versions = ch_versions.mix(PROCESS_FEATURE.out.versions)

        OPENMS_MZTABEXPORTER(PROCESS_FEATURE.out.consensusxml)
        ch_versions = ch_versions.mix(OPENMS_MZTABEXPORTER.out.versions)

    emit:
        consensusxml = PROCESS_FEATURE.out.consensusxml
        versions = ch_versions
}
