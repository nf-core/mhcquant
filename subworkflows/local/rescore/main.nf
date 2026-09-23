/*
 * Prepares the raw or compressed data holding spectra information for the subsequent database search.
 */

//
// MODULE: Loaded from modules/local/
//

include {
    MS2RESCORE ;
    MS2RESCORE as MS2RESCORE_GLOBAL
} from '../../../modules/local/ms2rescore'
include { OPENMS_PSMFEATUREEXTRACTOR                                  } from '../../../modules/nf-core/openms/psmfeatureextractor/main'
include {
    OPENMS_PERCOLATORADAPTER ;
    OPENMS_PERCOLATORADAPTER as OPENMS_PERCOLATORADAPTER_GLOBAL
} from '../../../modules/local/openmsthirdparty/percolatoradapter'
include {
    OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_GLOBAL ;
    OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_PSMS ;
    OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_PSMS_GLOBAL
} from '../../../modules/nf-core/openms/textexporter/main'
//
// MODULE: Installed directly from nf-core/modules
//

include { OPENMS_IDMERGER as OPENMS_IDMERGER_GLOBAL                   } from '../../../modules/nf-core/openms/idmerger/main'
include {
    OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE ;
    OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE_GLOBAL ;
    OPENMS_IDFILTER as OPENMS_IDFILTER_GLOBAL ;
    OPENMS_IDFILTER as OPENMS_IDFILTER_PSMS ;
    OPENMS_IDFILTER as OPENMS_IDFILTER_PSMS_GLOBAL
} from '../../../modules/nf-core/openms/idfilter/main'

workflow RESCORE {
    take:
    ch_merged_runs
    ch_multiqc_files

    main:
    // Compute features via ms2rescore. In ristretto mode the output already carries the
    // ristretto q-value as main score; in Percolator mode it carries the features for Percolator.
    MS2RESCORE(ch_merged_runs)

    if (params.rescoring_engine == 'ristretto') {
        ch_rescored_local = MS2RESCORE.out.idxml
        ch_global_input   = MS2RESCORE.out.idxml
    }
    else {
        // Read MS2Rescore feature names into meta so the nf-core module stays generic; -extra is set via ext.args
        MS2RESCORE.out.idxml
            .join(MS2RESCORE.out.feature_names)
            .map { meta, idxml, feature_names ->
                def extra = feature_names.readLines().drop(1).findAll { it.contains('\t') }.collect { it.split('\t', -1) }.findAll { !it[0].contains('psm_file') }.collect { it[1] }
                [meta + [extra_features: extra.join(' ')], idxml]
            }
            .set { ch_psmfeatureextractor_input }

        OPENMS_PSMFEATUREEXTRACTOR(ch_psmfeatureextractor_input)

        // Run Percolator with local FDR
        OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
        ch_multiqc_files = ch_multiqc_files.mix(OPENMS_PERCOLATORADAPTER.out.feature_weights.map { meta, feature_weights -> feature_weights })
        ch_rescored_local = OPENMS_PERCOLATORADAPTER.out.idxml
        ch_global_input   = OPENMS_PSMFEATUREEXTRACTOR.out.idxml.map { meta, idxml -> [meta.findAll { k, _v -> k != 'extra_features' }, idxml] }
    }

    if (params.global_fdr) {
        // Group by search_preset for global FDR. Samples without a preset all share
        // the same params (CLI or defaults), so they correctly group under 'global'.
        OPENMS_IDMERGER_GLOBAL(
            ch_global_input.map { group_meta, idxml -> [group_meta + [id: group_meta.search_preset ?: 'global'], idxml] }.groupTuple()
        )

        if (params.rescoring_engine == 'ristretto') {
            // Second MS²Rescore pass: all features are already present, so only ristretto runs, dataset-wide
            MS2RESCORE_GLOBAL(OPENMS_IDMERGER_GLOBAL.out.idxml.map { meta, idxml -> [meta, idxml, [], []] })
            ch_rescored_runs = MS2RESCORE_GLOBAL.out.idxml
        }
        else {
            // Run Percolator with global FDR (one per preset group)
            OPENMS_PERCOLATORADAPTER_GLOBAL(OPENMS_IDMERGER_GLOBAL.out.idxml)
            ch_rescored_runs = OPENMS_PERCOLATORADAPTER_GLOBAL.out.idxml
        }

        // Filter by global q-value
        OPENMS_IDFILTER_Q_VALUE_GLOBAL(ch_rescored_runs.map { id, idxml -> [id, idxml, []] })
        // Backfilter: match each local file with its corresponding preset's global FDR file
        OPENMS_IDFILTER_GLOBAL(
            ch_rescored_local.map { group_meta, idxml ->
                [group_meta.search_preset ?: 'global', group_meta, idxml]
            }.combine(
                OPENMS_IDFILTER_Q_VALUE_GLOBAL.out.filtered.map { global_meta, idxml -> [global_meta.id, idxml] },
                by: 0
            ).map { preset, group_meta, local_idxml, global_filtered_idxml ->
                [group_meta, local_idxml, global_filtered_idxml]
            }
        )
        ch_filter_q_value = OPENMS_IDFILTER_GLOBAL.out.filtered
        // Save globally merged runs in tsv (one per preset group)
        OPENMS_TEXTEXPORTER_GLOBAL(OPENMS_IDFILTER_Q_VALUE_GLOBAL.out.filtered)
        ch_global_rescored = ch_rescored_runs
        ch_global_filtered = OPENMS_IDFILTER_Q_VALUE_GLOBAL.out.filtered
    }
    else {
        ch_rescored_runs = ch_rescored_local
        // Filter by local q-value
        OPENMS_IDFILTER_Q_VALUE(ch_rescored_runs.map { group_meta, idxml -> [group_meta, idxml, []] })
        ch_filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.filtered
        ch_global_rescored = channel.empty()
        ch_global_filtered = channel.empty()
    }

    // PSM-level tables: every PSM of every FDR-passing peptidoform, taken from the pre-filter rescored
    // files (both engines keep all PSMs there) and whitelisted by the FDR-filtered output.
    OPENMS_IDFILTER_PSMS(ch_rescored_local.join(ch_filter_q_value))
    OPENMS_TEXTEXPORTER_PSMS(OPENMS_IDFILTER_PSMS.out.filtered)
    OPENMS_IDFILTER_PSMS_GLOBAL(ch_global_rescored.join(ch_global_filtered))
    OPENMS_TEXTEXPORTER_PSMS_GLOBAL(OPENMS_IDFILTER_PSMS_GLOBAL.out.filtered)

    ch_filter_q_value
        .map { meta, file -> [[id: meta.id], file] }
        .branch {
            // Empty FDR-filtered idXML (no peptides) is ~120 lines of OpenMS scaffolding.
            non_empty: it[1].countLines() > 130
            empty:     true
        }
        .set { ch_fdr_branched }

    emit:
    rescored_runs      = ch_rescored_runs.map { meta, file -> [[id: meta.id], file] }
    psms_tsv           = OPENMS_TEXTEXPORTER_PSMS.out.tsv.map { meta, file -> [[id: meta.id], file] }
    fdr_filtered       = ch_fdr_branched.non_empty
    fdr_filtered_empty = ch_fdr_branched.empty
    multiqc_files      = ch_multiqc_files
}
