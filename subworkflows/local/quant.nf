/*
 * Perform the quantification of the samples when the parameter --skip_quantification is not provided
 * This workflow splits the merged percolator output into the individual runs and filters them based on the q-value
 * It then aligns the retention times of the runs and merges the idxml files together to use them as id_ext in featurefinder
 * Finally, it performs the quantification and emits the consensusXML file
 */

include { OPENMS_IDRIPPER                          } from '../../modules/nf-core/openms/idripper/main'
include { OPENMS_IDSCORESWITCHER                   } from '../../modules/nf-core/openms/idscoreswitcher/main'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_QUANT } from '../../modules/nf-core/openms/idfilter/main'
include { OPENMS_IDMERGER as OPENMS_IDMERGER_QUANT } from '../../modules/nf-core/openms/idmerger/main'

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
        // Rip post-percolator idXML files and manipulate such that we end up with [meta_run1, idxml_run1] [meta_run2, idxml_run2] ...
        // We need to make sure that the order of the runs is the same as in the mzml files since IDRipper always sorts the runs
        // (and nextflow does not guarantee the order of the maps in merged_meta_map)
        OPENMS_IDRIPPER( merged_pout ).idxmls
                .flatMap { merged_meta, idxmls -> idxmls.collect { file -> [[spectra: file.baseName], file] } }
                // join on file basename to make sure that the order of the runs is the same as in the mzml files
                // Is there a smoother way to do this?
                .join( merge_meta_map
                        .flatMap { merged_meta, metas -> metas }
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
        ch_runs_score_switched
            .map { meta, idxml -> [[id: meta.sample + '_' + meta.condition], meta, idxml] }
            .combine(filter_q_value, by:0)
            .map { merge_id, meta, idxml, q_value -> [meta, idxml, q_value] }
            .set { ch_runs_to_filter}

        // Filter runs based on fdr filtered coprocessed percolator output.
        OPENMS_IDFILTER_QUANT( ch_runs_to_filter ).filtered
                .map { meta, idxml -> [[id:meta.sample + '_' + meta.condition], idxml] }
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
                                    .map { meta, aligned_idxml -> [[id: meta.sample + '_' + meta.condition], aligned_idxml] }
                                    .groupTuple())
        ch_versions = ch_versions.mix(OPENMS_IDMERGER_QUANT.out.versions)

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
        ch_versions = ch_versions.mix(PROCESS_FEATURE.out.versions)

    emit:
        consensusxml = PROCESS_FEATURE.out.consensusxml
        versions = ch_versions
}
