//
// Subworkflow with functionality specific to the nf-core/mhcquant pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { SDRF_CONVERT              } from '../../../modules/local/sdrf_convert/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input SDRF file

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Process SDRF file and convert to samplesheet format
    //
    SDRF_CONVERT(input)
    ch_versions = ch_versions.mix(SDRF_CONVERT.out.versions)

    //
    // Create channel from converted samplesheet
    //
    Channel
        .fromPath(SDRF_CONVERT.out.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def meta = [
                id: row.ID,
                sample: row.Sample,
                condition: row.Condition
            ]
            [ meta, file(row.ReplicateFileName) ]
        }
        .tap { ch_input }
        .groupTuple()
        // get number of files per sample-condition
        .map { meta, files -> [ meta, files.size()] }
        .combine( ch_input, by:0 )
        .map { meta, group_count, file -> [meta + ['group_count':group_count, 'spectra':file.baseName.tokenize('.')[0], 'ext':getCustomExtension(file)], file] }
        .set { ch_samplesheet }

    //
    // Create channel from the mandatory reference_database through params.fasta
    //
    Channel.fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [[id:fasta.getBaseName()], fasta] }
        .set { ch_fasta }

    emit:
    samplesheet = ch_samplesheet
    fasta       = ch_fasta
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Validate channels from input samplesheet
//
// Keeping this as an example for future samplesheet checks if additional fields are added (e.g. alleles)
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

def getCustomExtension(file) {
    def name = file.getName()
    if (name =~ /.*\.(d\.tar\.gz|d\.tar|d\.zip|mzML\.gz|raw|RAW|mzML|d)$/) {
        return name.split("\\.").drop(1).join(".").toLowerCase()
    } else {
        return file.getExtension().toLowerCase()
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    def citation_text = [
            "Tools used in the workflow included:",
            "OpenMS (Pfeuffer et al. 2024),",
            "DeepLC (Bouwmeester et al. 2021)",
            "MS²PIP (Declercq et al. 2023)",
            "MS²Rescore (Declercq et al. 2022)",
            "Percolator (Käll et al. 2007)",
            "MapAligner (Weisser et al. 2013)",
            "FeatureFinder (Weisser et al. 2017)",
            "MultiQC (Ewels et al. 2016)",
            "SDRF-Pipelines (BigBio)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    def reference_text = [
            "<li>Pfeuffer, J., Bielow, C., Wein, S. et al. OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data. Nat Methods (2024). doi: /10.1038/s41592-024-02197-7.</li>",
            "<li>Eng JK., Hoopman MR., Jahan, TA. et al. A Deeper Look into Comet—Implementation and Features. J. Am. Soc. Mass Spectrom. 2015, 26, 11, 1865–1874 (2015). doi: /10.1007/s13361-015-1179-x.</li>",
            "<li>Bouwmeester, R., Gabriels, R., Hulstaert, N. et al. DeepLC can predict retention times for peptides that carry as-yet unseen modifications. Nat Methods 18, 1363–1369 (2021). doi: /10.1038/s41592-021-01301-5<li>",
            "<li>Declercq A, Bouwmeester R, Chiva C, et al. Updated MS²PIP web server supports cutting-edge proteomics applications. Nucleic Acids Res. 2023 Jul 5;51(W1):W338-W342. doi: /10.1093/nar/gkad335<li>",
            "<li>Declercq A, Bouwmeester R, Hirschler A, Carapito C et al. MS2Rescore: Data-Driven Rescoring Dramatically Boosts Immunopeptide Identification Rates. Mol Cell Proteomics. 2022 Aug;21(8):100266. doi: /10.1016/j.mcpro.2022.100266<li>",
            "<li>Käll, L., Canterbury, J., Weston, J. et al. Semi-supervised learning for peptide identification from shotgun proteomics datasets. Nat Methods 4, 923–925 (2007). doi: /10.1038/nmeth1113<li>",
            "<li>Hendrik Weisser, Sven Nahnsen, Jonas Grossmann et al. An Automated Pipeline for High-Throughput Label-Free Quantitative Proteomics. Journal of Proteome Research 2013 12 (4), 1628-1644. doi: 10.1021/pr300992u<li>",
            "<li>Hendrik Weisser and Jyoti S. Choudhary, Journal of Proteome Research 2017 16 (8), 2964-2974. doi: /10.1021/acs.jproteome.7b00248<li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
            "<li>SDRF-Pipelines: A set of pipelines for extracting, processing and validating SDRF files. https://github.com/bigbio/sdrf-pipelines</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
