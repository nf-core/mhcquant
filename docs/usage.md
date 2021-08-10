# nf-core/mhcquant: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Samplesheet input](#samplesheet-input)
    * [Multiple runs of the same sample](#multiple-runs-of-the-same-sample)
    + [Full samplesheet](#full-samplesheet)
* [Running the pipeline](#running-the-pipeline)
    * [Updating the pipeline](#updating-the-pipeline)
    * [Reproducibility](#reproducibility)
* [Core Nextflow arguments](#core-nextflow-arguments)
    * [`-profile`](#-profile)
    * [`-resume`](#-resume)
* [Custom configuration](#custom-configuration)
* [Running in the background](#running-in-the-background)
* [Nextflow memory requirements](#nextflow-memory-requirements)

# nf-core/mhcquant: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mhcquant/usage](https://nf-co.re/mhcquant/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a tab-separated file with 4 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have multiple runs. The `Condition` identifiers can be used to further distinguish the sample groups.
Below is an example for the same sample sequenced across 3 lanes:

```console
ID\tSample\tCondition\tReplicateFileName
1\tWT\tA\t/path/to/MS/files/WT_A1.raw
2\tWT\tA\t/path/to/MS/files/WT_A2.raw
3\tWT\tA\t/path/to/MS/files/WT_A3.raw

```

### Full samplesheet

The pipeline will auto-detect whether a sample is either a mzML or raw files using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
ID\tSample\tCondition\tReplicateFileName
1\tWT\tA\t/path/to/MS/files/WT_A1.raw
2\tWT\tA\t/path/to/MS/files/WT_A2.raw
3\tWT\tA\t/path/to/MS/files/WT_A3.raw
4\tWT\tB\t/path/to/MS/files/WT_B1.raw
5\tWT\tB\t/path/to/MS/files/WT_B2.raw
6\tWT\tB\t/path/to/MS/files/WT_B3.raw
7\tKO\tA\t/path/to/MS/files/KO_A1.raw
8\tKO\tA\t/path/to/MS/files/KO_A2.raw
9\tKO\tA\t/path/to/MS/files/KO_A3.raw
10\tKO\tB\t/path/to/MS/files/KO_B1.raw
11\tKO\tB\t/path/to/MS/files/KO_B2.raw
12\tKO\tB\t/path/to/MS/files/KO_B3.raw

```

| Column         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `ID`       | An incrementing value which acts as a unique number for the given sample |
| `Sample`      | Custom sample name. This entry will be identical for multiple MS runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `Condition`      | Additional information of the sample can be defined here.|
| `ReplicateFileName` | Full path to the MS outcome file. These files have the extentions ".raw" or ".mzML" |

An [example samplesheet](../assets/samplesheet.tsv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:


```console
nextflow run nf-core/mhcquant --input 'samples.tsv' --fasta 'SWISSPROT_2020.fasta' --allele_sheet 'alleles.tsv' --vcf_sheet 'variants.tsv' --include_proteins_from_vcf --predict_class_1 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/mhcquant
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mhcquant releases page](https://github.com/nf-core/mhcquant/releases) and find the latest version number - numeric only (eg. `1.3`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/mhcquant`](https://hub.docker.com/r/nfcore/mhcquant/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/mhcquant`](https://hub.docker.com/r/nfcore/mhcquant/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

#### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

## Custom configuration

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
