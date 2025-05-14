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
        ch_versions = Channel.empty()

        // Compute group-wise alignment rt transformation
        OPENMS_MAPALIGNERIDENTIFICATION( ch_runs_to_be_aligned )
        ch_versions = ch_versions.mix(OPENMS_MAPALIGNERIDENTIFICATION.out.versions)

        // Join run specific trafoXMLs with meta information
        merge_meta_map
            .flatMap { group_meta, metas -> metas }
            .map { meta -> [[spectra:meta.spectra], meta]}
            .join( OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
                    .flatMap { group_meta, trafoxmls -> trafoxmls.collect { trafoxml -> [[spectra: trafoxml.baseName], trafoxml] } })
            .map { spectra, meta, trafoxml -> [meta, trafoxml] }
            .set { ch_trafos }

        // Align mzML files using trafoXMLs
        ch_trafos_mzmls = ch_mzml.join(ch_trafos)
        OPENMS_MAPRTTRANSFORMERMZML(ch_trafos_mzmls)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERMZML.out.versions)

        // Align idXMLfiles using trafoXMLs
        ch_runs_to_be_aligned
                .flatMap { group_meta, idxmls -> idxmls.collect { idxml -> [[spectra: idxml.baseName.replace("_fdr_filtered","")], idxml] } }
                .join( merge_meta_map
                        .flatMap { group_meta, metas -> metas }
                        .map { meta -> [[spectra:meta.spectra], meta]} )
                .map { group_meta, idxml, meta -> [meta, idxml] }
                .join( ch_trafos )
                .set { ch_trafos_idxml }

        OPENMS_MAPRTTRANSFORMERIDXML(ch_trafos_idxml)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERIDXML.out.versions)

    emit:
        versions = ch_versions
        aligned_idxml = OPENMS_MAPRTTRANSFORMERIDXML.out.aligned
        aligned_mzml = OPENMS_MAPRTTRANSFORMERMZML.out.aligned
}
