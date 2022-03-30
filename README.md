# ![nf-core/mhcquant](docs/images/nf-core-mhcquant_logo_light.png#gh-light-mode-only) ![nf-core/mhcquant](docs/images/nf-core-mhcquant_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/mhcquant/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/mhcquant/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/mhcquant/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/mhcquant/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/mhcquant/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.1569909-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.1569909)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23mhcquant-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/mhcquant)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nfcore/mhcquant** is a bioinformatics analysis pipeline used for quantitative processing of data dependent (DDA) peptidomics data.

It was specifically designed to analyse immunopeptidomics data, which deals with the analysis of affinity purified, unspecifically cleaved peptides that have recently been discussed intensively in [the context of cancer vaccines](https://www.nature.com/articles/ncomms13404).

The workflow is based on the OpenMS C++ framework for computational mass spectrometry. RAW files (mzML) serve as inputs and a database search (Comet) is performed based on a given input protein database. FDR rescoring is applied using Percolator based on a competitive target-decoy approach (reversed decoys). For label free quantification all input files undergo identification based retention time alignment (MapAlignerIdentification), and targeted feature extraction matching ids between runs (FeatureFinderIdentification). In addition, a variant calling file (vcf) can be specified to translate variants into proteins that will be included in the database search and binding predictions on specified alleles (alleles.tsv) using MHCFlurry (Class 1) or MHCNugget (Class 2) can be directly run on the output peptide lists. Moreover, if a vcf file was specified, neoepitopes will automatically be determined and binding predictions can also directly be predicted for them.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/mhcquant/results).

## Pipeline summary

1. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

![overview](assets/MHCquant_scheme.png)
(This chart was created with the help of [Lucidchart](https://www.lucidchart.com))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/mhcquant -profile test,YOURPROFILE
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run nf-core/mhcquant -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute> \
                                 --input 'samples.tsv' \
                                 --fasta 'SWISSPROT_2020.fasta' \
                                 --allele_sheet 'alleles.tsv' \
                                 --predict_class_1 \
                                 --refine_fdr_on_predicted_subset \
                                 --outdir ./results
   ```

## Documentation

The nf-core/mhcquant pipeline comes with documentation about the pipeline [usage](https://nf-co.re/mhcquant/usage), [parameters](https://nf-co.re/mhcquant/parameters) and [output](https://nf-co.re/mhcquant/output).

## Credits

nf-core/mhcquant was originally written by [Leon Bichmann](https://github.com/Leon-Bichmann) from the [Kohlbacher Lab](https://kohlbacherlab.org/). The pipeline was re-written in Nextflow DSL2 and is primarily maintained by [Marissa Dubbelaar](https://github.com/marissaDubbelaar) from [Clinical Collaboration Unit Translational Immunology](https://www.medizin.uni-tuebingen.de/en-de/das-klinikum/einrichtungen/kliniken/medizinische-klinik/kke-translationale-immunologie) and [Quantitative Biology Center](https://uni-tuebingen.de/forschung/forschungsinfrastruktur/zentrum-fuer-quantitative-biologie-qbic/) in Tübingen.

Helpful contributors:

- [Lukas Heumos](https://github.com/Zethson)
- [Alexander Peltzer](https://github.com/apeltzer)
- [Maxime Garcia](https://github.com/maxulysse)
- [Gisela Gabernet](https://github.com/ggabernet)
- [Susanne Jodoin](https://github.com/SusiJo)
- [Oskar Wacker](https://github.com/WackerO)
- [Leon Kuchenbecker](https://github.com/lkuchenb)
- [Phil Ewels](https://github.com/ewels)
- [Christian Fufezan](https://github.com/fu)
- [Sven Fillinger](https://github.com/sven1103)
- [Kevin Menden](https://github.com/KevinMenden)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#mhcquant` channel](https://nfcore.slack.com/channels/mhcquant) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/mhcquant for your analysis, please cite it using the following doi: [10.5281/zenodo.1569909](https://doi.org/10.5281/zenodo.1569909) and the corresponding manuscript:

> **MHCquant: Automated and Reproducible Data Analysis for Immunopeptidomics**
>
> Leon Bichmann, Annika Nelde, Michael Ghosh, Lukas Heumos, Christopher Mohr, Alexander Peltzer, Leon Kuchenbecker, Timo Sachsenberg, Juliane S. Walz, Stefan Stevanović, Hans-Georg Rammensee & Oliver Kohlbacher
>
> _Journal of Proteome Research_ 2019 18 (11), 3876-3884. doi: [10.1021/acs.jproteome.9b00313](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00313)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, references of tools and data used in this pipeline are as follows:

> **Fred2 Immunoinformatics Toolbox**
>
> Schubert B. et al, _Bioinformatics_ 2016 Jul 1;32(13):2044-6. doi: [10.1093/bioinformatics/btw113](https://academic.oup.com/bioinformatics/article/32/13/2044/1743767) Epub 2016 Feb 26
>
> **Comet Search Engine**
>
> Eng J.K. et al, _J Am Soc Mass Spectrom._ 2015 Nov;26(11):1865-74. doi: [10.1007/s13361-015-1179-x](https://pubs.acs.org/doi/10.1007/s13361-015-1179-x) Epub 2015 Jun 27.
>
> **Percolator**
>
> Käll L. et al, _Nat Methods_ 2007 Nov;4(11):923-5. doi: [10.1038/nmeth1113](https://www.nature.com/articles/nmeth1113) Epub 2007 Oct 21.
>
> **Identification based RT Alignment**
>
> Weisser H. et al, _J Proteome Res._ 2013 Apr 5;12(4):1628-44. doi: [10.1021/pr300992u](https://pubs.acs.org/doi/10.1021/pr300992u) Epub 2013 Feb 22.
>
> **Targeted peptide quantification**
>
> Weisser H. et al, _J Proteome Res._ 2017 Aug 4;16(8):2964-2974. doi: [10.1021/acs.jproteome.7b00248](https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00248) Epub 2017 Jul 19.
>
> **MHC affinity prediction**
>
> O'Donnell T.J., _Cell Syst._ 2018 Jul 25;7(1):129-132.e4. doi: [10.1016/j.cels.2018.05.014](https://www.sciencedirect.com/science/article/pii/S2405471218302321) Epub 2018 Jun 27.
>
> Shao X.M., _Cancer Immunol Res._ 2020 Mar;8(3):396-408. doi: [10.1158/2326-6066.CIR-19-0464](https://aacrjournals.org/cancerimmunolres/article/8/3/396/469597/High-Throughput-Prediction-of-MHC-Class-I-and-II) Epub 2019 Dec 23.
