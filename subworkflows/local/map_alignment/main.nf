/*
 * Align retention times of runs to be able to quantify them.
 */
include { OPENMS_MAPALIGNERIDENTIFICATION }                                 from '../../../modules/local/openms/mapaligneridentification'
include {
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERMZML
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERIDXML }               from '../../../modules/local/openms/maprttransformer'


workflow MAP_ALIGNMENT {
    take:
        ch_runs_to_be_aligned
        ch_mzml
        merge_meta_map

    main:
        // Compute group-wise alignment rt transformation
        OPENMS_MAPALIGNERIDENTIFICATION( ch_runs_to_be_aligned )

        // Run-specific trafoXMLs: [[spectra], trafoxml]
        OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
            .flatMap { group_meta, trafoxmls -> [trafoxmls].flatten().collect { trafoxml -> [[spectra: trafoxml.baseName], trafoxml] } }
            .set { ch_trafos }

        // Runs with their meta, idXML and mzML: [[spectra], meta, idxml, mzml]
        ch_runs_to_be_aligned
            .flatMap { group_meta, idxmls -> [idxmls].flatten().collect { idxml -> [[spectra: idxml.baseName.replace("_fdr_filtered","")], idxml] } }
            .join( merge_meta_map
                    .flatMap { group_meta, metas -> metas }
                    .map { meta -> [[spectra:meta.spectra], meta]} )
            .join( ch_mzml.map { meta, mzml -> [[spectra: meta.spectra], mzml] } )
            // Groups whose alignment failed (MapAlignerIdentification exit 8 is ignored) have no trafoXMLs
            .join( ch_trafos, remainder: true )
            .map { spectra, idxml, meta, mzml, trafoxml -> [meta, idxml, mzml, trafoxml] }
            .branch { meta, idxml, mzml, trafoxml ->
                aligned: trafoxml
                failed:  true
            }
            .set { ch_runs_by_alignment }

        // Groups that could not be aligned are excluded from quantification: [[id: sample_condition]]
        ch_runs_by_alignment.failed
            .map { meta, idxml, mzml, trafoxml -> [id: "${meta.sample}_${meta.condition}"] }
            .unique()
            .set { ch_failed_groups }

        ch_failed_groups.subscribe { group_meta ->
            log.warn "RT alignment of sample '${group_meta.id}' failed: at least one run has no peptide IDs in common with the other runs within --max_rt_alignment_shift. Skipping quantification for this sample, only identifications are reported."
        }

        // Align mzML and idXML files using trafoXMLs
        OPENMS_MAPRTTRANSFORMERMZML(ch_runs_by_alignment.aligned.map { meta, idxml, mzml, trafoxml -> [meta, mzml, trafoxml] })
        OPENMS_MAPRTTRANSFORMERIDXML(ch_runs_by_alignment.aligned.map { meta, idxml, mzml, trafoxml -> [meta, idxml, trafoxml] })

    emit:
        aligned_idxml = OPENMS_MAPRTTRANSFORMERIDXML.out.aligned
        aligned_mzml  = OPENMS_MAPRTTRANSFORMERMZML.out.aligned
        trafoxml      = OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
        failed_groups = ch_failed_groups
}
