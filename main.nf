#!/usr/bin/env nextflow

/*
========================================================================================
                            nf-core/mhcquant
========================================================================================
    nf-core/mhcquant Analysis Pipeline.
    #### Homepage / Documentation
    https://github.com/nf-core/mhcquant
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    log.info nfcoreHeader()
    def command = "nextflow run nf-core/mhcquant --input 'sample_sheet.tsv' --fasta 'SWISSPROT_2020.fasta'  --allele_sheet 'alleles.tsv'  --predict_class_1  --refine_fdr_on_predicted_subset -profile standard,docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
        include { MHCQUANT } from './workflows/mhcquant' addParams( summary_params: summary_params )
        MHCQUANT ()
}
