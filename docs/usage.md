# nf-core/mhcquant: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/mhcquant/usage](https://nf-co.re/mhcquant/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## SDRF input

You will need to create an SDRF (Sample and Data Relationship Format) file with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. SDRF is a standardized format used in proteomics that contains rich metadata about samples.

```bash
--input '[path to SDRF file]'
```

### SDRF format

SDRF files contain detailed sample metadata in a tabular format with specific columns for sample characteristics, data file locations, and experimental factors. The pipeline uses sdrf-pipelines to process this file and extract the necessary information.

Key columns in an SDRF file include:

| Column                      | Description                                                                                           |
| --------------------------- | ----------------------------------------------------------------------------------------------------- |
| `source name`               | Custom sample name. This entry will be identical for multiple MS runs from the same sample.           |
| `characteristics[organism]` | The organism from which the sample was derived                                                        |
| `comment[data file]`        | Full path to the MS file. These files have the extentions .raw, .mzML, mzML.gz, .d, .d.tar.gz, .d.zip |
| `factor value[...]`         | Experimental factors that can be used for grouping samples                                            |
| `characteristics[mhc type]` | MHC alleles present in the sample (important for MHC peptide analysis)                                |

The pipeline will auto-detect whether a sample is either in mzML, raw or tdf file format using the information provided in the SDRF file.

An [example SDRF file](../assets/example_sdrf.tsv) has been provided with the pipeline.

### Multiple runs of the same sample

MS runs are merged based on the `source name` and factor values before they are rescored with `Percolator`. Typically technical replicates of a sample are merged together to report one peptide list per sample.

For example, in the provided example SDRF file, there are multiple technical replicates for each sample (e.g., sample1_classI, sample1_classII), and these will be merged together during processing.

### Creating SDRF files

You can create SDRF files manually or use tools like the [SDRF Composer](https://github.com/bigbio/sdrf-composer) to generate them. The sdrf-pipelines tool also provides validation functionality to ensure your SDRF file is correctly formatted:

```bash
parse_sdrf.py validate -s your_sdrf_file.tsv
```

## Recommended search settings

Fine-tuning search settings is important to obtain the most optimal results for your MS data. _These settings heavily depend on the MS instrument settings used to generate the data_. If you want to reprocess public data, make sure you use the settings mentioned in the methods section! The following table acts as an orientation of commonly used search settings for instruments:

| MS-Device                | timsTOF  |          | Orbitrap Fusion Lumos |          | Q Exactive Orbitrap |          | LTQ Orbitrap XL |          |
| ------------------------ | -------- | -------- | --------------------- | -------- | ------------------- | -------- | --------------- | -------- |
|                          | class I  | class II | class I               | class II | class I             | class II | class I         | class II |
| instrument               | high_res | high_res | high_res              | high_res | high_res            | high_res | low_res         | low_res  |
| digest_mass_range        | 800:2500 | 800:5000 | 800:2500              | 800:5000 | 800:2500            | 800:5000 | 800:2500        | 800:5000 |
| activation_method        | CID      | CID      | HCD                   | HCD      | HCD                 | HCD      | CID             | CID      |
| prec_charge              | 1:4      | 1:5      | 2:3                   | 2:5      | 2:3                 | 2:5      | 2:3             | 2:5      |
| precursor_error_units    | ppm      | ppm      | ppm                   | ppm      | ppm                 | ppm      | ppm             | ppm      |
| number_mods              | 3        | 5        | 3                     | 5        | 3                   | 5        | 3               | 5        |
| precursor_mass_tolerance | 20       | 20       | 5                     | 5        | 5                   | 5        | 5               | 5        |
| fragment_mass_tolerance  | 0.02     | 0.02     | 0.02                  | 0.02     | 0.02                | 0.02     | 0.50025         | 0.50025  |
| fragment_bin_offset      | 0        | 0        | 0                     | 0        | 0                   | 0        | 0.4             | 0.4      |

Modifications are specified via `--variable_mods` and `fixed_mods` using the [UNIMOD nomenclature](https://www.unimod.org/unimod_help.html) via OpenMS. Check out [helper page](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_CometAdapter.html) of OpenMS for the full list of options. Multiple modifications are specified as `'Oxidation (M),Acetyl (N-term),Phospho (S)'`.

Further information about the command line arguments is documented on the [nf-core website](https://nf-co.re/mhcquant/dev/parameters) or by using `--help`.

## Rescoring using MS²Rescore

By default the pipline generates additional features using MS²PIP and DeepLC via the MS²Rescore framework (`--feature_generators deeplc,ms2pip`). Additional feature generators can be added (`basic,deeplc,ms2pip,ionmob,im2deep`) to boost identification rates and quality. Please make sure you provide the correct `--ms2pip_model` (default: `Immuno-HCD`). All available MS²PIP models can be found on [GitHub](https://github.com/compomics/ms2pip).

MS²Rescore creates a comprehensive QC report of the added features used for rescoring. This report is only available if `--rescoring_engine mokapot` is specified (default: `percolator`). The report can be found in `<OUTDIR>/multiqc/ms2rescore`. Further information on the tool itself can be read up in the published paper [Declerq et al. 2022](<https://www.mcponline.org/article/S1535-9476(22)00074-3/fulltext>)

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/mhcquant \
  --input 'sdrf_file.tsv' \
  --outdir <OUTDIR> \
  --fasta 'SWISSPROT_2020.fasta' \
  <SEARCH PARAMS> \
  --peptide_min_length 8 \
  --peptide_max_length 14 \
  --ms2pip_model 'Immuno-HCD' \
  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/mhcquant -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './sdrf_file.tsv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/mhcquant
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/mhcquant releases page](https://github.com/nf-core/mhcquant/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

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

```bash
NXF_OPTS='-Xms1g -Xmx4g'
