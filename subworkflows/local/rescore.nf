/*
 * Prepares the raw or compressed data holding spectra information for the subsequent database search.
 */

//
// MODULE: Loaded from modules/local/
//

include { MS2RESCORE                                                  } from '../../modules/local/ms2rescore'
include { OPENMS_PSMFEATUREEXTRACTOR                                  } from '../../modules/local/openms_psmfeatureextractor'
include { OPENMS_PERCOLATORADAPTER;
          OPENMS_PERCOLATORADAPTER as OPENMS_PERCOLATORADAPTER_GLOBAL } from '../../modules/local/openms_percolatoradapter'
include { OPENMS_TEXTEXPORTER as OPENMS_TEXTEXPORTER_GLOBAL           } from '../../modules/local/openms_textexporter'
//
// MODULE: Installed directly from nf-core/modules
//

include { OPENMS_IDMERGER as OPENMS_IDMERGER_GLOBAL                   } from '../../modules/nf-core/openms/idmerger/main'
include { OPENMS_IDSCORESWITCHER                                      } from '../../modules/nf-core/openms/idscoreswitcher/main.nf'
include { OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE;
          OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE_GLOBAL;
          OPENMS_IDFILTER as OPENMS_IDFILTER_GLOBAL                   } from '../../modules/nf-core/openms/idfilter/main'

workflow RESCORE {

    take:
        ch_merged_runs

    main:
        ch_versions = Channel.empty()

    // Compute features via ms2rescore
    MS2RESCORE(ch_merged_runs)
    ch_versions = ch_versions.mix(MS2RESCORE.out.versions)

    if (params.rescoring_engine == 'mokapot') {
        log.warn "The rescoring engine is set to mokapot. This rescoring engine currently only supports psm-level-fdr via ms2rescore."
        // Switch comet e-value to mokapot q-value
        OPENMS_IDSCORESWITCHER(MS2RESCORE.out.idxml)
        ch_versions = ch_versions.mix(OPENMS_IDSCORESWITCHER.out.versions)
        ch_rescored_runs = OPENMS_IDSCORESWITCHER.out.idxml

    } else {
		// Extract PSM features for Percolator
		OPENMS_PSMFEATUREEXTRACTOR(MS2RESCORE.out.idxml.join(MS2RESCORE.out.feature_names))
		ch_versions = ch_versions.mix(OPENMS_PSMFEATUREEXTRACTOR.out.versions)

        // Run Percolator with local FDR
		OPENMS_PERCOLATORADAPTER(OPENMS_PSMFEATUREEXTRACTOR.out.idxml)
		ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions)
        ch_pout = OPENMS_PERCOLATORADAPTER.out.idxml

		if (params.global_fdr) {
            // Merge all samples into one group
			OPENMS_IDMERGER_GLOBAL(OPENMS_PSMFEATUREEXTRACTOR.out.idxml.map {group_meta, idxml -> [[id:'global'], idxml] }.groupTuple())
            // Run Percolator with global FDR
            OPENMS_PERCOLATORADAPTER_GLOBAL(OPENMS_IDMERGER_GLOBAL.out.idxml)
            ch_rescored_runs = OPENMS_PERCOLATORADAPTER_GLOBAL.out.idxml
            // Filter by global percolator q-value
            OPENMS_IDFILTER_Q_VALUE_GLOBAL(ch_rescored_runs.map {id, idxml -> [id, idxml, []]})
            ch_versions = ch_versions.mix(OPENMS_IDFILTER_Q_VALUE_GLOBAL.out.versions)
            // Backfilter sample_condition runs according to global FDR
            OPENMS_IDFILTER_GLOBAL(ch_pout.combine(OPENMS_IDFILTER_Q_VALUE_GLOBAL.out.filtered.map{ it[1] }))
            ch_filter_q_value = OPENMS_IDFILTER_GLOBAL.out.filtered
            // Save globally merged runs in tsv
            OPENMS_TEXTEXPORTER_GLOBAL(OPENMS_PERCOLATORADAPTER_GLOBAL.out.idxml)

		} else {
            ch_rescored_runs = ch_pout
            // Filter by percolator q-value
            OPENMS_IDFILTER_Q_VALUE(ch_rescored_runs.map {group_meta, idxml -> [group_meta, idxml, []]})
            ch_versions = ch_versions.mix(OPENMS_IDFILTER_Q_VALUE.out.versions)
            ch_filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.filtered
        }
	}

    emit:
        rescored_runs = ch_rescored_runs
        fdr_filtered = ch_filter_q_value
        versions = ch_versions
}
