/*
 * Perform the quantification of the samples when the parameter --skip_quantification is not provided
 */
include { OPENMS_MAPALIGNERIDENTIFICATION }                                 from '../../modules/local/openms_mapaligneridentification'
include {
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERMZML
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERIDXML }               from '../../modules/local/openms_maprttransformer'


workflow MAP_ALIGNMENT {
    take:
        fdr_filtered_runs
        mzml_files

    main:
        ch_versions = Channel.empty()
        // Compute group-wise alignment rt transformation
        OPENMS_MAPALIGNERIDENTIFICATION(fdr_filtered_runs.groupTuple())
        ch_versions = ch_versions.mix(OPENMS_MAPALIGNERIDENTIFICATION.out.versions.first().ifEmpty(null))
        // MapAligner module might return trafoXML files and meta information in the wrong order, so we need to take care of that
        joined_trafos = OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
            .transpose()
            .flatMap {
                meta, trafoxml ->
                    ident = trafoxml.baseName.split('_').first()
                    [[[id:ident, sample:meta.sample, condition:meta.condition, ext:meta.ext], trafoxml]]
            }
        // Intermediate step to join RT transformation files with mzml channels -> [meta, idxml, mzml]
        joined_trafos_mzmls = mzml_files.join(joined_trafos)
        // Intermediate step to join RT transformation files with idxml channels -> [meta, idxml, trafoxml]
        joined_trafos_ids = fdr_filtered_runs
                                    .map {
                                        merge_id, meta, idxml -> [meta, idxml]
                                    }
                                    .join(joined_trafos)
        // Align mzML files using trafoXMLs
        OPENMS_MAPRTTRANSFORMERMZML(joined_trafos_mzmls)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERMZML.out.versions.first().ifEmpty(null))
        // Align unfiltered idXMLfiles using trafoXMLs
        OPENMS_MAPRTTRANSFORMERIDXML(joined_trafos_ids)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERIDXML.out.versions.first().ifEmpty(null))

    emit:
        // Define the information that is returned by this workflow
        versions = ch_versions
        aligned_idxml = OPENMS_MAPRTTRANSFORMERIDXML.out.aligned
        aligned_mzml = OPENMS_MAPRTTRANSFORMERMZML.out.aligned
}
