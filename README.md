# ![nf-core/mhcquant](docs/images/nf-core-mhcquant_logo.png)

**Identify and quantify peptides from mass spectrometry raw data**.

[![GitHub Actions CI Status](https://github.com/nf-core/mhcquant/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/mhcquant/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/mhcquant/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/mhcquant/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/mhcquant.svg)](https://hub.docker.com/r/nfcore/mhcquant)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23mhcquant-4A154B?logo=slack)](https://nfcore.slack.com/channels/mhcquant)

## Introduction

nfcore/mhcquant is a bioinformatics analysis pipeline used for quantitative processing of data dependent (DDA) peptidomics data.

It was specifically designed to analyse immunopeptidomics data, which deals with the analysis of affinity purified, unspecifically cleaved peptides that have recently been discussed intensively in [the context of cancer vaccines](https://www.nature.com/articles/ncomms13404).

The workflow is based on the OpenMS C++ framework for computational mass spectrometry. RAW files (mzML) serve as inputs and a database search (Comet) is performed based on a given input protein database. FDR rescoring is applied using Percolator based on a competitive target-decoy approach (reversed decoys). For label free quantification all input files undergo identification based retention time alignment (MapAlignerIdentification), and targeted feature extraction matching ids between runs (FeatureFinderIdentification). In addition, a variant calling file (vcf) can be specified to translate variants into proteins that will be included in the database search and binding predictions on specified alleles (alleles.tsv) using MHCFlurry (Class 1) or MHCNugget (Class 2) can be directly run on the output peptide lists. Moreover, if a vcf file was specified, neoepitopes will automatically be determined and binding predictions can also directly be predicted for them.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

![overview](assets/MHCquant_scheme.png)
(This chart was created with the help of [Lucidchart](https://www.lucidchart.com))

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/mhcquant -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/mhcquant -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
                                  --input 'samples.tsv'
                                  --fasta 'SWISSPROT_2020.fasta'
                                  --allele_sheet 'alleles.tsv' 
                                  --predict_class_1 
                                  --refine_fdr_on_predicted_subset
    ```

See [usage docs](https://nf-co.re/mhcquant/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:

* Protein database addition by mutated genome variants (`Fred2 Immunoinformatics Toolbox`)
* Database search ('Comet')
* False discovery rate estimation ('Percolator')
* Retention time alignment ('OpenMS-MapAlignerIdentification')
* Targeted peptide quantification ('OpenMS-FeatureFinderIdentification')
* MHC peptide affinity prediction ('MHCFlurry','MHCNuggets')

## Documentation

The nf-core/mhcquant pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/mhcquant/usage) and [output](https://nf-co.re/mhcquant/output).

## Credits

We thank the following people for their extensive assistance in the development
of this pipeline:

* Leon Bichmann (@Leon-Bichmann)
* Lukas Heumos (@Zethson)
* Alexander Peltzer (@apeltzer)

The pipeline was converted to Nextflow DSL2 by Marissa Dubbelaar (@marissaDubbelaar)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#mhcquant` channel](https://nfcore.slack.com/channels/mhcquant) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use `nf-core/mhcquant` for your analysis, please cite:

> **MHCquant: Automated and Reproducible Data Analysis for Immunopeptidomics**
>
> Leon Bichmann, Annika Nelde, Michael Ghosh, Lukas Heumos, Christopher Mohr, Alexander Peltzer, Leon Kuchenbecker, Timo Sachsenberg, Juliane S. Walz, Stefan Stevanović, Hans-Georg Rammensee & Oliver Kohlbacher
>
> Journal of Proteome Research 2019 18 (11), 3876-3884
> DOI: 10.1021/acs.jproteome.9b00313

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, references of tools and data used in this pipeline are as follows:

> **Fred2 Immunoinformatics Toolbox**
>
> Schubert B. et al, _Bioinformatics_ 2016 Jul 1;32(13):2044-6. doi: 10.1093/bioinformatics/btw113. Epub 2016 Feb 26
>
> **Comet Search Engine**
>
> Eng J.K. et al, _J Am Soc Mass Spectrom._ 2015 Nov;26(11):1865-74. doi: 10.1007/s13361-015-1179-x. Epub 2015 Jun 27.
>
> **Percolator**
>
> Käll L. et al, _Nat Methods_ 2007 Nov;4(11):923-5. doi: 10.1038/nmeth1113. Epub 2007 Oct 21.
>
> **Identification based RT Alignment**
>
> Weisser H. et al, _J Proteome Res._ 2013 Apr 5;12(4):1628-44. doi: 10.1021/pr300992u. Epub 2013 Feb 22.
>
> **Targeted peptide quantification**
>
> Weisser H. et al, _J Proteome Res._ 2017 Aug 4;16(8):2964-2974. doi: 10.1021/acs.jproteome.7b00248. Epub 2017 Jul 19.
>
> **MHC affinity prediction**
>
> O'Donnell T.J., _Cell Syst._ 2018 Jul 25;7(1):129-132.e4. doi: 10.1016/j.cels.2018.05.014. Epub 2018 Jun 27.
>
> Shao X.M., _Cancer Immunol Res._ 2020 Mar;8(3):396-408. doi: 10.1158/2326-6066.CIR-19-0464. Epub 2019 Dec 23.
