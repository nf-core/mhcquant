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
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { SDRF_TO_SAMPLESHEET       } from '../sdrf_to_samplesheet'

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
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:


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

    def before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/mhcquant ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    def after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/mhcquant/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Detect input type and build samplesheet channels
    //
    def inputType = detectInputType(params.input)

    if (params.qpx_out && (!params.quantify || !(inputType in ['sdrf', 'pride_id']))) {
        error("--qpx_out requires both --quantify and an SDRF/PRIDE input (--input <SDRF file> or PXD accession).")
    }

    if (inputType == 'sdrf' || inputType == 'pride_id') {
        //
        // SDRF / PRIDE input: samplesheet is produced by a process, so validate lazily.
        //
        def sdrf_path  = (inputType == 'sdrf') ? params.input : null
        def pride_id   = (inputType == 'pride_id') ? params.input : null

        if (inputType == 'sdrf' && !pride_id) {
            def matcher = (file(params.input).name =~ /PXD\d{6,}/)
            pride_id = matcher.find() ? matcher.group() : null
        }

        SDRF_TO_SAMPLESHEET(sdrf_path, pride_id)

        ch_presets_file = SDRF_TO_SAMPLESHEET.out.search_presets
        ch_samplesheet_rows = SDRF_TO_SAMPLESHEET.out.samplesheet
            .flatMap { samplesheet_path ->
                samplesheetToList(samplesheet_path.toString(), "${projectDir}/assets/schema_input.json")
            }

    } else {
        //
        // Standard samplesheet: parse eagerly so validation fails fast with a documented error.
        //
        ch_presets_file = channel.fromPath(params.search_presets, checkIfExists: true)
        ch_samplesheet_rows = channel.fromList(
            samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
        )
    }

    //
    // Build presets map (shared)
    //
    ch_presets_map = ch_presets_file
        .map { presets_file ->
            samplesheetToList(presets_file.toString(), "${projectDir}/assets/schema_search_presets.json")
                .collectEntries { item ->
                    // samplesheetToList wraps all-meta rows in a list
                    def row = (item instanceof List) ? item[0] : item
                    // nf-schema parses empty TSV cells as [] instead of ''; normalize for string operations
                    ['fixed_mods', 'variable_mods'].each { key ->
                        def v = row[key]
                        if (!v || (v instanceof String && !v.trim())) row[key] = ''
                    }
                    [(row.preset_name): row]
                }
        }

    //
    // Enrich rows and resolve search params (shared)
    //
    ch_samplesheet_rows
        .map { meta, file, fasta ->
            def m = meta + [sample: meta.sample.toString(), condition: meta.condition.toString()]
            [m.subMap('sample', 'condition'), m, file, fasta]
        }
        .tap { ch_input }
        .groupTuple()
        .map { group_meta, metas, files, fastas -> [group_meta, files.size()] }
        .combine(ch_input, by: 0)
        .map { group_meta, group_count, meta, file, fasta ->
            def enrichedMeta = meta + [group_count: group_count, spectra: file.baseName.tokenize('.')[0], ext: getCustomExtension(file)]
            [enrichedMeta, file, fasta]
        }
        .set { ch_samplesheet_raw }

    ch_samplesheet = ch_samplesheet_raw
        .combine(ch_presets_map)
        .map { meta, file, fasta, presetsMap ->
            [resolvePresetParams(meta, presetsMap), file]
        }

    //
    // Create channel from the reference_database through params.fasta
    //
    if (params.fasta) {
        channel.fromPath(params.fasta, checkIfExists: true)
            .map { fasta -> [[id:fasta.getBaseName()], fasta] }
            .set { ch_fasta }

        ch_samplesheet_raw
            .map{ meta, file, fasta -> fasta }
            .flatten()
            .first()
            .subscribe {
                log.warn """\
                    Both --fasta and samplesheet FASTA files were provided!
                    The pipeline will use --fasta (${params.fasta}), ignoring samplesheet FASTA entries.
                    To use the samplesheet FASTA files instead, remove the --fasta parameter.
                    """.stripIndent()
            }
    } else {
        // Fasta from samplesheet column
        ch_fasta = ch_samplesheet_raw.map { meta, file, fasta -> [groupKey([id: "${meta.sample}_${meta.condition}"], meta.group_count), fasta] }
        ch_fasta
            .map { meta, fasta -> fasta }
            .flatten()
            .ifEmpty {
                error '''\
                    Error: No FASTA files provided.
                    Please either:
                    1. Use --fasta parameter, or
                    2. Include a 'Fasta' column in your samplesheet
                    '''.stripIndent()
            }
        ch_fasta
            .groupTuple()
            .map { group_meta, fastas -> [group_meta, fastas.first()] }
            .set { ch_fasta }
    }

    emit:
    samplesheet = ch_samplesheet
    fasta       = ch_fasta
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

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Detect input type from the --input parameter value
//
def detectInputType(input) {
    def inputStr = input.toString()
    if (inputStr =~ /^PXD\d{6,}$/) {
        return 'pride_id'
    }
    if (inputStr.endsWith('.sdrf.tsv')) {
        return 'sdrf'
    }
    return 'samplesheet'
}

//
// Resolve search parameters for a row. When the row sets a preset, the preset's values
// win for every key the preset defines (the preset is sealed and cannot be overridden
// via --<param>, -params-file, or -c). Keys not defined by the preset (or rows without
// a preset) fall through to params[key], which Nextflow resolves through its native
// precedence (CLI > -params-file > config defaults).
//
def resolvePresetParams(meta, presetsMap) {
    def searchParamKeys = ['instrument', 'activation_method', 'digest_mass_range', 'prec_charge',
                           'precursor_mass_tolerance', 'precursor_error_units', 'fragment_mass_tolerance',
                           'fragment_bin_offset', 'number_mods', 'ms2pip_model',
                           'peptide_min_length', 'peptide_max_length',
                           'fixed_mods', 'variable_mods']
    def presetName = meta.search_preset
    def hasPreset = presetName && !(presetName instanceof List && presetName.size() == 0) && presetName != ''
    def presetConfig = hasPreset ? presetsMap[presetName] : [:]
    if (hasPreset && !presetConfig) {
        error "Unknown search preset '${presetName}'. Available: ${presetsMap.keySet().join(', ')}"
    }
    if (!presetConfig) { presetConfig = [:] }

    def result = new LinkedHashMap(meta)
    searchParamKeys.each { key ->
        result.put(key, presetConfig.containsKey(key) ? presetConfig[key] : params[key])
    }
    return result
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
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
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
