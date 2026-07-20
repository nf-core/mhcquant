/*
 * Perform the quantification of the samples when the parameter --quantify is provided
 * This workflow splits the merged percolator output into the individual runs and filters them based on the q-value
 * It then aligns the retention times of the runs and merges the idxml files together to use them as id_ext in featurefinder
 * Finally, it performs the quantification and emits the consensusXML file
 */

include { OPENMS_IDRIPPER                          } from '../../../modules/nf-core/openms/idripper/main'
include { OPENMS_IDSCORESWITCHER                   } from '../../../modules/nf-core/openms/idscoreswitcher/main'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_QUANT } from '../../../modules/nf-core/openms/idfilter/main'
include { OPENMS_IDMERGER as OPENMS_IDMERGER_QUANT } from '../../../modules/nf-core/openms/idmerger/main'
include { OPENMS_MZTABEXPORTER                     } from '../../../modules/local/openms/mztabexporter'

include { MAP_ALIGNMENT                            } from '../map_alignment'
include { PROCESS_FEATURE                          } from '../process_feature'

workflow QUANT {
    take:
        merge_meta_map
        merged_pout
        filter_q_value
        mzml

    main:
        // Split post-percolator idXML files and manipulate such that we end up with [meta_run1, idxml_run1] [meta_run2, idxml_run2] ...
        // We need to make sure that the order of the runs is the same as in the mzml files since IDRipper always sorts the runs
        // (and nextflow does not guarantee the order of the maps in merged_meta_map)
        OPENMS_IDRIPPER( merged_pout ).idxmls
            // Handle both single file and list of files
            .flatMap { group_meta, idxmls -> [idxmls].flatten().collect { idxml -> [[spectra: idxml.baseName], idxml] } }
            // join on file basename to make sure that the order of the runs is the same as in the mzml files
            .join( merge_meta_map
                    .flatMap { group_meta, metas -> metas }
                    .map { meta -> [[spectra:meta.spectra], meta]} )
            .map { spectra, idxmls, meta -> [meta, idxmls] }
            .set { ch_ripped_idxml }

        // Switch to xcorr for filtering since q-values are set to 1 with peptide-level-fdr
        if (params.fdr_level == 'peptide_level_fdrs'){
            ch_runs_score_switched = OPENMS_IDSCORESWITCHER( ch_ripped_idxml ).idxml
        } else {
            ch_runs_score_switched = ch_ripped_idxml
        }

        // Manipulate such that [meta_run1, idxml_run1, pout_group1], [meta_run2, idxml_run2, pout_group1] ...
        ch_runs_score_switched
            .map { meta, idxml -> [groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count) , meta, idxml] }
            // Unwrap to plain [id] map: asymmetric GroupKey.equals() makes combine() drop pairs by arrival order (nextflow-io/nextflow#4104)
            .map { key, meta, idxml -> [key.getGroupTarget(), meta, idxml] }
            .combine(filter_q_value, by:0)
            .map { group_meta, meta, idxml, q_value -> [meta, idxml, q_value] }
            .set { ch_runs_to_filter}

        // Filter runs based on fdr filtered coprocessed percolator output.
        OPENMS_IDFILTER_QUANT( ch_runs_to_filter ).filtered
                .map { meta, idxml -> [ groupKey([id:"${meta.sample}_${meta.condition}"], meta.group_count), idxml] }
                .groupTuple()
                .set { ch_runs_to_be_aligned }

        // Align retention times of runs
        MAP_ALIGNMENT(
            ch_runs_to_be_aligned,
            mzml,
            merge_meta_map
        )

        // We need to merge groupwise the aligned idxml files together to use them as id_ext in featurefinder
        OPENMS_IDMERGER_QUANT( MAP_ALIGNMENT.out.aligned_idxml
                                    .map { meta, aligned_idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), aligned_idxml] }
                                    .groupTuple())

        // Manipulate channels such that we end up with : [meta, mzml, run_idxml, merged_runs_idxml]
        MAP_ALIGNMENT.out.aligned_mzml
                .join( MAP_ALIGNMENT.out.aligned_idxml )
                .map { meta, aligned_mzml, idxml -> [ groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), meta, aligned_mzml, idxml] }
                .groupTuple()
                .map { group_meta, meta, aligned_mzml, idxml -> [group_meta, meta, aligned_mzml, idxml] }
                .join( OPENMS_IDMERGER_QUANT.out.idxml )
                .map { group_meta, meta, aligned_mzml, idxml, merged_idxml -> [meta, aligned_mzml, idxml, merged_idxml] }
                .transpose()
                .set { ch_runs_to_be_quantified }

        PROCESS_FEATURE ( ch_runs_to_be_quantified )

        OPENMS_MZTABEXPORTER(PROCESS_FEATURE.out.consensusxml)

    emit:
        consensusxml = PROCESS_FEATURE.out.consensusxml
        trafoxml = MAP_ALIGNMENT.out.trafoxml
}
