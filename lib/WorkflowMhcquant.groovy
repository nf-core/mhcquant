//
// This file holds several functions specific to the workflow/rnaseq.nf in the nf-core/rnaseq pipeline
//

class WorkflowMhcquant {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        // Error raise when the database fasta is not included
        if (!params.fasta) { 
            log.error "No database fasta was provided, make sure you have used the '--fasta' option."
            System.exit(1)
        }

        // Warning for the outpout dir (when not defined)
        if (!params.outdir) {
            outdirWarn(log)
        }

        // Database options
        if (params.skip_decoy_generation) { 
            skipDecoyWarn(log)
        }

        // FDR Scoring
        if (params.klammer && params.description_correct_features == 0) {
            klammerConflictWarn(log)
        }

        // Quantification options
        if (params.quantification_fdr) {
            enabledFDRWarn(log)
        }
    }


    //
    // Print a warning when the output directory is undefined
    //
    private static void outdirWarn(log) {
        log.warn "=============================================================================\n" +
            "  Results into \'./results\'.\n" +
            "  If you want to define a result directory, please use the --outdir option" +
            "============================================================================="
    }

    //
    // Print a warning when the parameter skip_decoy_generation is used
    //
    private static void skipDecoyWarn(log) {
        log.warn "=============================================================================\n" +
            "  Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement\n" +
            "  Decoys have to be named with DECOY_ as prefix in your fasta database\n" +
            "==================================================================================="
    }

    //
    // Print a warning when the parameter skip_decoy_generation is used
    //
    private static void klammerConflictWarn(log) {
        log.warn "=============================================================================\n" +
            "  Quantification FDR enabled\n" + 
            "==================================================================================="
    }

    //
    // Print a warning when the parameter skip_decoy_generation is used
    //
    private static void enabledFDRWarn(log) {
        log.warn "=============================================================================\n" +
            "  Klammer was specified, but description of correct features was still 0. Please provide a description of correct features greater than 0.\n" + 
            "  Klammer has been turned off!\n" + 
            "==================================================================================="
    }

}
