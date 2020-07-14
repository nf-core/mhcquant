# nf-core/mhcquant: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the final output of the pipeline should include the following files:

* [all_features_merged_resolved.mzTab](#mzTab) - the community standard format for sharing mass spectrometry search results
* [all_features_merged_resolved.csv](#csv) - aggregate csv report, containing all information about peptide identification and quantification results
* [found_neoepitopes.csv](#found_neoepitopes) - a csv listing all neoepitopes found in the mass spectrometry search, independant of binding predictions
* [vcf_neoepitopes.csv](#vcf_neoepitopes) - a csv listing all theoretically possible neoepitope sequences from the variants specified in the vcf
* [_vcf.fasta](#fasta) - the fasta database including mutated proteins used for the database search
* [class_1/2_binding_predictions](#class_1/2_binding_predictions) - a folder containing the respective binding predictions of all detected peptides and all theoretically possible neoepitope sequences
* [Intermediate_resuls](#intermediates) - a folder containing all intermediate results from the steps in the pipeline (unfiltered and filtered PSMs, aligned mzMLs, features, etc. ..)
* [Intermediate_results/mhcquant_file_order.txt](#file_order) - a txt file listing the order of files used for annotating intensitites in the mzTab output
* [Documentation](#docs) - a folder containing summarized reports of the pipeline execution
* [pipeline_info](#info) - a folder containing detailed reports on computational runtimes and workflow steps

## mzTab

The output mzTab contains many columns annotating the most important information - here are a few outpointed:

```bash
PEP   sequence   accession   best_search_engine_score[1]   retention_time   charge   mass_to_charge   peptide_abundance_study_variable[1]
```

Most important to know that in this format we annotated the q-value of each peptide identification in the best_seach_engine_score[1] column and peptide quantities in the peptide_abundance_study_variable columns.

[mzTab](http://www.psidev.info/mztab) is a light-weight format to report mass spectrometry search results. It provides all important information about identified peptide hits and is compatible with the PRIDE Archive - proteomics data repository:

Griss, J. et al. The mzTab Data Exchange Format: Communicating Mass-spectrometry-based Proteomics and Metabolomics Experimental Results to a Wider Audience. Mol Cell Proteomics 13, 2765-2775 (2014)

## csv

The csv output file is a table containing all information extracted from a database search throughout the pipeline. See the OpenMS or PSI documentation for more information about [annotated scores and format](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_TextExporter.html).

```bash
#MAP    id      filename        label   size
```

MAP contains information about the different mzML files that were provided initially

```bash
#RUN    run_id  score_type      score_direction date_time       search_engine_version   parameters
```

RUN contains information about the search that was performed on each run

```bash
#PROTEIN        score   rank    accession       protein_description     coverage        sequence
```

PROTEIN contains infomration about the protein ids corresponding to the peptides that were detected (No protein inference was performed)

```bash
#UNASSIGNEDPEPTIDE      rt      mz      score   rank    sequence        charge  aa_before       aa_after        score_type      search_identifier       accessions      FFId_category   feature_id      file_origin     map_index       spectrum_reference      COMET:IonFrac   COMET:deltCn    COMET:deltLCn   COMET:lnExpect  COMET:lnNumSP   COMET:lnRankSP  MS:1001491      MS:1001492      MS:1001493      MS:1002252      MS:1002253      MS:1002254      MS:1002255      MS:1002256      MS:1002257      MS:1002258      MS:1002259      num_matched_peptides    protein_references      target_decoy
```

UNASSIGNEDPEPTIDE contains information about PSMs that were identified but couldn't be quantified to a precursor feature on MS Level 1.

```bash
#CONSENSUS      rt_cf   mz_cf   intensity_cf    charge_cf       width_cf        quality_cf      rt_0    mz_0    intensity_0     charge_0        width_0 rt_1    mz_1    intensity_1     charge_1        width_1 rt_2    mz_2    intensity_2     charge_2        width_2 rt_3    mz_3    intensity_3     charge_3        width_3
```

CONSENSUS contains information about precursor features that were identified in multiple runs (eg. run 1-3 in this case)

```bash
#PEPTIDE        rt      mz      score   rank    sequence        charge  aa_before       aa_after        score_type      search_identifier       accessions      FFId_category   fea
```

PEPTIDE contains information about peptide hits that were identified and correspond to the consensus features described one row above.

## found_neoepitopes

csv file listing detected neoepitope sequences:

```bash
peptide sequence   geneID
```

## vcf_neoepitopes

csv file listing theoretically possible neoepitope sequences:

```bash
Sequence        Antigen ID       Variants
```

## class_1/2_binding_predictions

The prediction outputs are comma separated table (csv) for each allele, listing each peptide sequence and its corresponding predicted affinity scores:

```bash
peptide   allele   prediction   prediction_low   prediction_high   prediction_percentile
```
