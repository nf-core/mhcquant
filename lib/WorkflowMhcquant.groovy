//
// This file holds several functions specific to the workflow/mhcquant.nf in the nf-core/mhcquant pipeline
//

import groovy.text.SimpleTemplateEngine

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

        // Allele sheet
        if ( ( params.predict_class_1 || params.predict_class_2 ) && !params.allele_sheet)  {
            // If no allele sheet is provided
            noAllelesError(log)
            System.exit(1)
        }

        if ( params.include_proteins_from_vcf && !params.vcf_sheet)  {
            // If no vcf sheet is provided
            noVcfError(log)
            System.exit(1)
        }
    }

    //
    // Print a warning when the output directory is undefined
    //
    private static void outdirWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Results into \'./results\'.\n" +
            "  If you want to define a result directory, please use the --outdir option" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning when the parameter skip_decoy_generation is used
    //
    private static void skipDecoyWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement\n" +
            "  Decoys have to be named with DECOY_ as prefix in your fasta database\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning when the parameter skip_decoy_generation is used
    //
    private static void klammerConflictWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Quantification FDR enabled\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print a warning when the parameter skip_decoy_generation is used
    //
    private static void enabledFDRWarn(log) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Klammer was specified, but description of correct features was still 0\n."+
            "  Please provide a description of correct features greater than 0.\n" +
            "  Klammer has been turned off!\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print an error when the no allele_sheet is provided
    //
    private static void noAllelesError(log) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  --predict_class_1 or --predict_class_2 was used but no allele sheet has been provided. \n" +
            "  Make sure you have used the '--allele_sheet' option and include an allele sheet.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }

    //
    // Print an error when there is no vcf_sheet provided
    //
    private static void noVcfError(log) {
        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  --include_proteins_from_vcf was used but no vcf sheet has been provided. \n" +
            "  Make sure you have used the '--vcf_sheet' option and include an vcf sheet.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
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

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }
}
