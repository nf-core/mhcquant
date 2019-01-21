# nf-core/mhcquant: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
    * [`--mzmls`](#--mzmls)
    * [`--fasta`](#--fasta)
    * [`-profile`](#-profile-single-dash)
        * [`docker`](#docker)
        * [`awsbatch`](#awsbatch)
        * [`standard`](#standard)
        * [`none`](#none)
* [Mass Spectrometry Search](#Mass Spectrometry Search)
    * [`--peptide_min_length`](#--peptide_min_length)
    * [`--peptide_max_length`](#--peptide_max_length)
    * [`--fragment_mass_tolerance`](#--fragment_mass_tolerance)
    * [`--precursor_mass_tolerance`](#--precursor_mass_tolerance)
    * [`--fragment_bin_offset`](#--fragment_bin_offset)
    * [`--fdr_threshold`](#--fdr_threshold)
    * [`--fdr_level`](#--fdr_level)
    * [`--number_mods`](#--number_mods)
    * [`--num_hits`](#--num_hits)
    * [`--digest_mass_range`](#--digest_mass_range)
    * [`--pick_ms_levels`](#--pick_ms_levels)
    * [`--centroided`](#--centroided)
    * [`--prec_charge`](#--prec_charge)
    * [`--digest_mass_range`](#--digest_mass_range)
    * [`--activation_method`](#--activation_method)
    * [`--enzyme`](#--enzyme)
    * [`--fixed_mods`](#--fixed_mods)
    * [`--variable_mods`](#--variable_mods)
    * [`--spectrum_batch_size`](#--spectrum_batch_size)
* [Optional binding predicion](#optional-binding-prediction)
    * [`--run_prediction`](#--run_prediction)
    * [`--refine_fdr_on_predicted_subset`](#--refine_fdr_on_predicted_subset)
    * [`--affinity_threshold_subset`](#--affinity_threshold_subset)
    * [`--alleles`](#--alleles)
* [Optional variant translation](#optional-variant_translation)
    * [`--include_proteins_from_vcf`](#--include_proteins_from_vcf)
    * [`--vcf`](#--vcf)
    * [`--variant_annotation_style`](#--variant_annotation_style)
    * [`--variant_reference`](#--variant_reference)
    * [`--variant_indel_filter`](#--variant_indel_filter)
    * [`--variant_frameshift_filter`](#--variant_frameshift_filter)
    * [`--variant_snp_filter`](#--variant_snp_filter)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
    * [`-awsbatch`](#-awsbatch)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_email`](#--plaintext_email)


## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run nf-core/mhcquant --mzmls '*.mzML' --fasta 'SWISSPROT_12_2018.fasta' --alleles 'alleles.tsv' --vcf 'variants.vcf' -profile standard,docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/mhcquant
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mhcquant releases page](https://github.com/nf-core/mhcquant/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main Arguments

### `--mzmls`
Use this to specify the location of your input mzML files. For example:

```bash
--mzmls 'path/to/data/*.mzML'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character

### `--fasta`
If you prefer, you can specify the full path to your fasta input protein database when you run the pipeline:

```bash
--fasta '[path to Fasta protein database]'
```

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile standard,docker` - the order of arguments is important!

* `standard`
    * The default profile, used if `-profile` is not specified at all.
    * Runs locally and expects all software to be installed and available on the `PATH`.
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`nfcore/mhcquant`](http://hub.docker.com/r/nfcore/mhcquant/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `binac`
    * A profile for the BinAC cluster at the University of Tübingen 
    * Pulls software from Docker Hub via Singularity
* `cfc`
    * A profile for the Core Facility Cluster (CFC) at QBiC Tübingen
    * Pulls software from Docker Hub via Singularity
* `awsbatch`
    * A generic configuration profile to be used with AWS Batch.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

## Mass Spectrometry Search
### `--peptide_min_length`
Specify the minimum length of peptides considered after processing

### `--peptide_max_length`
Specify the maximum length of peptides considered after processing

### `--fragment_mass_tolerance`
Specify the fragment mass tolerance used for the comet database search. For High-Resolution instruments a fragment mass tolerance value of 0.02 is recommended. (See the Comet parameter documentation: eg. 0.02)

### `--precursor_mass_tolerance`
Specify the precursor mass tolerance used for the comet database search. For High-Resolution instruments a precursor mass tolerance value of 5ppm is recommended. (eg. 5)

### `--fragment_bin_offset`
Specify the fragment bin offset used for the comet database search. For High-Resolution instruments a fragment bin offset of 0 is recommended. (See the Comet parameter documentation: eg. 0)

### `--fdr_threshold`
Specify the false discovery rate threshold at which peptide hits should be selected. (eg. 0.01)

### `--fdr_level`
Specify the level at which the false discovery rate should be computed. 'peptide-level-fdrs' is recommended. ('peptide-level-fdrs', 'psm-level-fdrs', 'protein-level-fdrs')

### `--number_mods`
Specify the maximum number of modifications that should be contained in a peptide sequence match. (eg. 3)

### `--num_hits`
Specify the number of hits that should be reported for each spectrum. (eg. 1)

### `--digest_mass_range`
Specify the mass range that peptides should fullfill to be considered for peptide spectrum matching. (eg. 800:2500)

### `--pick_ms_levels`
If one ms level in the raw ms data is not centroided, specify the level here. (eg. 2)

### `--centroided`
Choose whether the specified ms_level in pick_ms_levels is centroided or not. ("True", "False")

### `--prec_charge`
Specifiy the precursor charge range that peptides should fullfill to be considered for peptide spectrum matching. (eg. "2:3")

### `--activation method`
Specify which fragmentation method was used in the MS acquisition ('ALL', 'CID', 'ECD', 'ETD', 'PQD', 'HCD', 'IRMPD')

### `--enzyme`
Specify which enzymatic restriction should be applied ('unspecific cleavage', 'Trypsin', see OpenMS enzymes)

### `--fixed_mods`
Specify which fixed modifications should be applied to the database search (eg. '' or 'Carbamidomethyl (C)', see OpenMS modifications)

### `--variable_mods`
Specify which variable modifications should be applied to the database search (eg. 'Oxidation (M)', see OpenMS modifications)

### `--spectrum_batch_size`
Size of Spectrum batch for Comet processing (Decrease/Increase depending on Memory Availability)

## Optional binding prediction
### `--run_prediction`
Set to 'True' or 'False' depending on whether binding predictions using the tool mhcflurry should be run. (Check whether your alleles are supported by mhcflurry)

### `--refine_fdr_on_predicted_subset`
Set to 'True' or 'False' depending on whether binding predictions using the tool mhcflurry should be used to subset all PSMs not passing the q-value threshold. If specified the FDR will be refined using Percolator on the subset of predicted binders among all PSMs resulting in an increased identification rate. (Check whether your alleles are supported by mhcflurry)

### `--affinity_threshold_subset`
Affinity threshold (nM) used to define binders for PSM subset selection in the fdr refinement procedure (eg. 500)

### `--alleles`
Specify a .tsv file containing the alleles of your probes. (line separated)

## Optional variant translation
### `--include_proteins_from_vcf`
Set to 'True' or 'False' depending on whether variants should be translated to proteins and included into your fasta for database search.

### `--vcf`
Specify a .vcf file containing the information about genomic variants.

### `--variant_annotation_style`
Specify style of tool used for variant annotation - currently supported: "SNPEFF", "VEP", "ANNOVAR"

### `--variant_reference`
Specify genomic reference used for variant annotation: "GRCH37", "GRCH38"

### `--variant_indel_filter`
Specify whether insertions and deletions should not be considered for variant translation

### `--variant_frameshift_filter`
Specify whether frameshifts should not be considered for variant translation

### `--variant_snp_filter`
Specify whether snps should not be considered for variant translation

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, you can specify a config file using `-c` that contains the following:

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

