# nf-core/mhcquant: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## General

<summary>Output files</summary>

- `*.tsv`
- `*.mzTab` (if `--quantify` is specified)
- `global_fdr/global.tsv` (if `--global_fdr` is specified)
- `spectrum_library/*_speclib.tsv` (if `--generate_speclib` is specified)

#### TSV
The TSV output of MHCquant is a tab-delimited file holding information about FDR-filtered peptides and their properties such as retention time, charge, mass-to-charge, and protein accessions. Additionally, a selected number of features produced by `Comet` and `MS²Rescore` are stored. Understanding these scores and features is crucial in narrowing down peptides of interest.

| **Tool**   | **Metric**                      | **Description**                                                                                   |
|------------|---------------------------------|---------------------------------------------------------------------------------------------------|
| Comet      | `COMET:xcorr`                   | Cross-correlation score between observed and theoretical spectra; higher scores indicate better matches. |
| Comet      | `COMET:deltaCn`                 | Normalized difference between the top two XCorr scores; higher values indicate better discrimination between top matches. |
| Comet      | `COMET:deltaLCn`                | Similar to deltaCn but considers the difference between the top match and the best match with a different peptide sequence. |
| Comet      | `COMET:lnExpect`                | Natural logarithm of the Expectation value; lower values suggest more statistically significant matches. |
| DeepLC     | `predicted_retention_time_best` | Predicted retention time for the precursor peak (all PSMs with the same sequence, modifications, and charge) closest to the predicted retention time. |
| DeepLC     | `observed_retention_time_best`  | Observed retention time for the precursor peak (all PSMs with the same sequence, modifications, and charge) closest to the predicted retention time. |
| DeepLC     | `rt_diff_best`                  | Difference between observed and predicted retention time for the precursor peak (all PSMs with the same sequence, modifications, and charge) closest to the predicted retention time. |
| MS2PIP     | `spec_pearson`                  | Pearson correlation coefficient between observed and predicted b and y ion intensities of MS2 spectra; higher values indicate better agreement. |
| MS2PIP     | `std_abs_diff`                  | Standard deviation of absolute differences between observed and predicted b and y ion intensities; lower values indicate better agreement. |
| IM2Deep    | `ccs_observed_im2deep`          | Observed peptide collisional cross-section (CCS) value. |
| IM2Deep    | `ccs_predicted_im2deep`         | Predicted peptide CCS value by IM2Deep. |
| IM2Deep    | `ccs_error_im2deep`             | Difference between observed and predicted CCS values; indicates prediction accuracy. |

The TSV file in quantification mode (by using `--quantify`) additionally holds the columns `rt_cf|0|1|..`, `mz_cf|0|1|..`, `intensity_cf|0|1|..`, and `charge_cf|0|1|..`, where `cf` (consensus feature) represents the aggregated value of the individual MS runs (`0|1|..`) that are part of the `<Sample>_<Condition>` group in the input samplesheet. Use e.g. the individual intensities of MS runs to compare two groups and conduct a differential-presentation analysis.

#### mzTab
The mzTab output file follows the a [HUPO-PSI format](<https://www.mcponline.org/article/S1535-9476(20)32821-8/fulltext>) and combines all information of the sample-condition group extracted from a database search throughout the pipeline. A detailed explanation of the respective entries are elaborately explained [here](https://psidev.info/sites/default/files/2017-07/R2_The_ten_minute_guide_to_mzTab.pdf). MzTab files are compatible with the PRIDE Archive - proteomics data repository and can be uploaded as search files.

MzTab files contain many columns and annotate the most important information - here are a few outpointed:

```bash
PEP  sequence  accession  best_search_engine_score[1]  retention_time  charge  mass_to_charge  peptide_abundance_study_variable[1]
```

By default (only identification) the `best_search_engine_score[1]` holds the percolator q-value. If `--quantify` is specified the Comet XCorr of each peptide identification is annotated in the `best_search_engine_score[1]` column and peptide quantities in the `peptide_abundance_study_variable` columns.


#### Global FDR
Global FDR controls the accumulation of false positives by estimating the proportion of incorrect identifications across the whole dataset, rather than per-spectrum or local identifications. This approach ensures that as more peptide-spectrum matches (PSMs) are accepted, the overall rate of erroneous entries remains statistically bounded, preventing inflation of false discoveries. If `--global_fdr` is specified, global FDR is applied on each group of `{Sample}_{Condition}` of the samplesheet. Additionally, a merged TSV of all MSruns with global FDR filter can then be found in the `global_fdr` folder.

#### Spectrum library
Experimental spectrum libraries from DDA experiments capture observed peptide fragmentation patterns and retention times. EasyPQP converts these into spectral libraries compatible with search engines like DIA-NN. Here is how a spectrum library looks like:
```tsv title="speclib.tsv
PrecursorMz	ProductMz	Annotation	ProteinId	GeneName	PeptideSequence	ModifiedPeptideSequence	PrecursorCharge	LibraryIntensity	NormalizedRetentionTime	PrecursorIonMobility
271.824095	118.086256	y1^1	sp|P62917|RL8_HUMAN	-	DAPAGRKV	DAPAGRKV	3	122.085594	18.952473326867114	1.12423
271.824095	187.071335	b2^1	sp|P62917|RL8_HUMAN	-	DAPAGRKV	DAPAGRKV	3	546.7115	18.952473326867114	1.12423
```

</details>

### Intermediate results

<details  markdown="1">

This folder contains the intermediate results from various steps of the MHCquant pipeline (e.g. (un)filtered PSMs, aligned mzMLs, features)

<summary>Output files</summary>
    
- `intermediate_results/`
  - `alignment`: Contains the `trafoXML` files of each run that document the retention time shift after alignment in quantification mode.
  - `comet`: Contains pin files generated by comet after database search
  - `features`: Holds information of quantified features in `featureXML` files as a result of the [FeatureFinderIdentification](https://openms.de/doxygen/release/3.0.0/html/TOPP_FeatureFinderIdentification.html) in the quantification mode.
  - `ion_annotations`(if `--annotate_ions` is specified)
    - `{Sample}_{Condition}_all_peaks.tsv`: Contains metadata of all measured ions of peptides reported after peptide identification.
    - `{Sample}_{Condition}_matching_ions.tsv`: Contains ion annotations and additional metadata of peptides reported after peptide identification.
  - `rescoring`
    - `{Sample}_{Condition}_(psm|ms2rescore).idXML`: File holding extra features generated by MS²Rescore that will be used by percolator or mokapot.
    - `{Sample}_{Condition}_pout.idXML`: Unfiltered percolator output.
    - `{Sample}_{Condition}_pout_filtered.idXML`: FDR-filtered percolator output.


</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`

- `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.

- `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.

- `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`

  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.html`.
  - Reports generated by the pipeline: `software_versions.yml`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
