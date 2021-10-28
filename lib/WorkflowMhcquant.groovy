//
// This file holds several functions specific to the workflow/mhcquant.nf in the nf-core/mhcquant pipeline
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

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

}
