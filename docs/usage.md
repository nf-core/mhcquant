# nf-core/mhcquant: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`--mzmls`](#--mzmls)
  * [`--raw_input`](#--raw_input)
  * [`--raw_files`](#--raw_files)
  * [`--fasta`](#--fasta)
  * [`-profile`](#-profile)
* [Mass Spectrometry Search](#Mass-Spectrometry-Search)
  * [`--peptide_min_length`](#--peptide_min_length)
  * [`--peptide_max_length`](#--peptide_max_length)
  * [`--fragment_mass_tolerance`](#--fragment_mass_tolerance)
  * [`--precursor_mass_tolerance`](#--precursor_mass_tolerance)
  * [`--fragment_bin_offset`](#--fragment_bin_offset)
  * [`--use_a_ions`](#--use_a_ions)
  * [`--use_c_ions`](#--use_c_ions)
  * [`--use_x_ions`](#--use_x_ions)
  * [`--use_z_ions`](#--use_z_ions)
  * [`--fdr_threshold`](#--fdr_threshold)
  * [`--fdr_level`](#--fdr_level)
  * [`--number_mods`](#--number_mods)
  * [`--num_hits`](#--num_hits)
  * [`--digest_mass_range`](#--digest_mass_range)
  * [`--pick_ms_levels`](#--pick_ms_levels)
  * [`--run_centroidisation`](#--run_centroidisation)
  * [`--prec_charge`](#--prec_charge)
  * [`--digest_mass_range`](#--digest_mass_range)
  * [`--activation_method`](#--activation_method)
  * [`--enzyme`](#--enzyme)
  * [`--fixed_mods`](#--fixed_mods)
  * [`--variable_mods`](#--variable_mods)
  * [`--max_rt_alignment_shift`](#--max_rt_alignment_shift)
  * [`--spectrum_batch_size`](#--spectrum_batch_size)
  * [`--skip_decoy_generation`](#--skip_decoy_generation)
  * [`--quantification_fdr`](#--quantification_fdr)
  * [`--quantification_min_prob`](#--quantification_min_prob)
  * [`--predict_RT`](#--predict_RT)
* [Optional binding predicion](#optional-binding-prediction)
  * [`--predict_class_1`](#--predict_class_1)
  * [`--predict_class_2`](#--predict_class_2)
  * [`--refine_fdr_on_predicted_subset`](#--refine_fdr_on_predicted_subset)
  * [`--affinity_threshold_subset`](#--affinity_threshold_subset)
  * [`--class_1_alleles`](#--class_1_alleles)
  * [`--class_2_alleles`](#--class_2_alleles)
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
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/mhcquant --mzmls '*.mzML' --fasta 'SWISSPROT_12_2018.fasta' --class_1_alleles 'alleles.tsv' --vcf 'variants.vcf' --include_proteins_from_vcf --predict_class_1 -profile docker
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

First, go to the [nf-core/mhcquant releases page](https://github.com/nf-core/mhcquant/releases) and find the latest version number - numeric only (eg. `1.3`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `--mzmls`

Use this to specify the location of your input mzML files. For example:

```bash
--mzmls 'path/to/data/*.mzML'
```

### `--raw_input`

Set this flag if you want to use raw files instead of mzml files

### `--raw_files`

Use this to specify the location of your input raw files. For example:

```bash
--raw_input --raw_files 'path/to/data/*.raw'
```

### `--fasta`

If you prefer, you can specify the full path to your fasta input protein database when you run the pipeline:

```bash
--fasta '[path to Fasta protein database]'
```

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/mhcquant`](http://hub.docker.com/r/nfcore/mhcquant/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/mhcquant`](http://hub.docker.com/r/nfcore/mhcquant/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

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

### `--use_a_ions`

Include a ions into the peptide spectrum matching

### `--use_c_ions`

Include c ions into the peptide spectrum matching

### `--use_x_ions`

Include x ions into the peptide spectrum matching

### `--use_z_ions`

Include z ions into the peptide spectrum matching

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

### `--run_centroidisation`

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

Multiple fixed or variable modifications can be specified comma separated (e.g. 'Carbamidomethyl (C),Oxidation (M)')

### `--max_rt_alignment_shift`

Set a maximum retention time shift for the linear rt alignment

### `--spectrum_batch_size`

Size of Spectrum batch for Comet processing (Decrease/Increase depending on Memory Availability)

### `--skip_decoy_generation`

If you want to use your own decoys, you can specify a databaset that includes decoy sequences. However, each database entry should keep the prefix 'DECOY_'.
One should consider though that this option will then prevent to append variants to the database and if not using reversed decoys the subset refinement FDR option will not work.

### `--quantification_fdr`

Set this option to assess and assign quantification of peptides with an FDR measure (Weisser H. and Choudhary J.S. J Proteome Res. 2017 Aug 4)

### `--quantification_min_prob`

Specify a cut off probability value for quantification events as a filter

### `--predict_RT`

Set this option to predict times of all identified peptides and possible neoepitopes based on high scoring ids

## Optional binding prediction

### `--predict_class_1`

Set flag depending on whether MHC class 1 binding predictions using the tool mhcflurry should be run. [Check whether your alleles are supported by mhcflurry](supported_class_1_alleles.md)

### `--predict_class_2`

Set flag depending on whether MHC class 2 binding predictions using the tool mhcnugget should be run. [Check whether your alleles are supported by mhcnugget](supported_class_2_alleles.md)

### `--refine_fdr_on_predicted_subset`

Set to 'True' or 'False' depending on whether binding predictions using the tool mhcflurry should be used to subset all PSMs not passing the q-value threshold. If specified the FDR will be refined using Percolator on the subset of predicted binders among all PSMs resulting in an increased identification rate. (Please be aware that this option is only available for MHC class I data of alleles that are supported by mhcflurry)

### `--affinity_threshold_subset`

Affinity threshold (nM) used to define binders for PSM subset selection in the fdr refinement procedure (eg. 500)

### `--class_1_alleles`

Specify a .tsv file containing the MHC class 1 alleles of your probes. (line separated)

### `--class_2_alleles`

Specify a .tsv file containing the MHC class 2 alleles of your probes. (line separated)

## Optional variant translation

### `--include_proteins_from_vcf`

Set to 'True' or 'False' depending on whether variants should be translated to proteins and included into your fasta for database search.

### `--vcf`

Specify a .vcf file containing the information about genomic variants (vcf < v.4.2).

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

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

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

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.
