# nf-core/mhcquant: Changelog

## v2.3.1 nfcore/mhcquant "White Gold Swallow" - 2022/05/10

### `Added`

- [#210](https://github.com/nf-core/mhcquant/pull/210) - Icons to the different parameters in the `nextflow_schema.json`

### `Fixed`

- [#211](https://github.com/nf-core/mhcquant/issues/211) Resolved the problem that there was no output from `OPENMS_MZTABEXPORTER_QUANT`
- [#212](https://github.com/nf-core/mhcquant/issues/212) - Altered the CometAdapter functionalities for resolve the issue with the `fixed_mods`

### `Dependencies`

### `Deprecated`

## v2.3.0 nfcore/mhcquant "White Gold Swallow" - 2022/04/05

### `Added`

- [#206](https://github.com/nf-core/mhcquant/issues/206) Updated the workflow picture
- Adjustments of the `PRE_QUANTIFICATION` subworkflow: `OPENMS_FALSEDISCOVERYRATE`, `OPENMS_IDFILTER_FOR_ALIGNMENT`, and `OPENMS_TEXTEXPORTER_PSMS`
- Included `OPENMS_TEXTEXPORTER_UNQUANTIFIED`to write a combined FDR filtered output file for unquantified data
- Included `pipeline summary` and increment the `documentation` paragraph
- [#195](https://github.com/nf-core/mhcquant/issues/195) Updated parameter documentation
- [#189](https://github.com/nf-core/mhcquant/issues/189) Added backslashes in Quick Start in README
- [#188](https://github.com/nf-core/mhcquant/issues/188) Added reference links to README

### `Fixed`

- Typo in previous release date
- [#208](https://github.com/nf-core/mhcquant/pull/208) - nf-core template update (version 2.3.2)
- [#199](https://github.com/nf-core/mhcquant/issues/199) Fixes some typos and stuff in the output documentation
- [#192](https://github.com/nf-core/mhcquant/issues/192) Fixed samplesheet format in usage.md
- [#184](https://github.com/nf-core/mhcquant/issues/184) Fix parsing for VEP annotated VCF files

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `MultiQC`  | 1.11        | 1.12        |

### `Deprecated`

- [#191](https://github.com/nf-core/mhcquant/issues/191) Removed the table of contents from usage.md

## v2.2.0 nfcore/mhcquant "Silver Titanium Ostrich" - 2022/01/14

### `Added`

- Included the newest nf-core template (version 2.2)
- Adjustment of the README, including all contributors
- Inclusion of the PSMs files (tsv format) per replicates in `results/PSMs`
- Include check in WorkflowMhcquant, to determine if the allele and vcf sheet has been provided under specific circumstances

### `Fixed`

- Changed parameters in the nextflow_schema.json to be in coherence with the nextflow.config
- Error that was raised in generate_proteins_from_vcf
- Problems that were detected in predict_possible_class1_neoepitopes and predict_possible_class2_neoepitopes
- Error that occurred in mhcnuggets_predictneoepitopesclass2 (faulty container set up)

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `fred2`    | 2.0.6       | 2.0.7       |

### `Deprecated`

## v2.1.0 nf-core/mhcquant "Olive Tin Hamster" - 2021/12/09

### `Added`

- Inclusion of assets/schema_input.json
- Added the multiQC again to report the versions
- MHCquant parameters are now directly assigned to the argument of the process

### `Fixed`

- Fixed typos
- [#165] - Raise memory requirements of FeatureFinderIdentification step
- [#176] - Pipeline crashes when setting the --skip_quantification flag

### `Dependencies`

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency            | Old version | New version |
| --------------------- | ----------- | ----------- |
| `openms`              | 2.5.0       | 2.6.0       |
| `openms-thirdparty`   | 2.5.0       | 2.6.0       |
| `thermorawfileparser` | 1.2.3       | 1.3.4       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

### `Deprecated`

## v2.0.0 nf-core/mhcquant "Steel Beagle" - 2021/09/03

### `Added`

DSL1 to DSL2 conversion

    - Different processes based on a unique step in the pipeline
    - Inclusion of one sub workflow: refine fdr on predicted subset
    - The process: openms_cometadapter includes commented lines (which could be used as a reference for future module development)
    - MHCquant pipeline is ran from workflows/mhcquant.nf instead of main
    - Template update to nf-core tools version 2.1

### `Fixed`

### `Dependencies`

### `Deprecated`

## v1.6.0 nf-core/mhcquant "Beijing Duck" - 2020/09/11

### `Added`

- Template raise to 1.10.2
- Added parameter json schema
- Added full size AWS test profile
- Included new parameters for Neutral loss and precursor ion inclusion

### `Fixed`

- Changed trigger for AWS tests

### `Dependencies`

### `Deprecated`

## v1.5.1 nf-core/mhcquant "Flying Fish" - 2020/04/24

### `Added`

### `Fixed`

- set optimal config for cluster execution
- fix duplication of ids / mixing of channels

### `Dependencies`

### `Deprecated`

## v1.5.0 nf-core/mhcquant "Flying Fish" - 2020/04/18

### `Added`

- integrate sample, allele and vcf sheets instead of file dirs
- branched mzML/raw input
- introduce param to skip quantification

### `Fixed`

- raise OpenMS version to 2.5
- adapt workflow accoringly with new options
- remove specifying input as file dirs eg "data/\*.mzML"

### `Dependencies`

### `Deprecated`

## v1.4.0 nf-core/mhcquant "Blue Kingfisher" - 2020/03/18

### `Added`

- Raw File Reading
- RT prediction
- Quantification FDR
- Variant pass filter
- nf-core template update 1.8 and 1.9
- Added version numbers of mhcnuggets and Fred2

### `Fixed`

- output file order in intermediate results
- increased run times for MS search and variant translation

### `Dependencies`

### `Deprecated`

## v1.3.0 nf-core/mhcquant "Red Parrot" - 2019/08/03

### `Added`

- nf-core template update
- x,z,a,c ions
- quantification fdri

### `Fixed`

- empty neoepitope list bugs fixed
- documentation
- scrape version numbers

### `Dependencies`

### `Deprecated`

## v1.2.6 nf-core/mhcquant "Golden Eagle" - 2019/05/05

### `Added`

- MHCnugget predictor
- Few fixes
- RT features for percolator

### `Fixed`

### `Dependencies`

### `Deprecated`

## v1.2.6 nf-core/mhcquant "Golden Eagle" - 2019/03/05

### `Added`

### `Fixed`

- linear retention time alignment
- refine_fdr README

### `Dependencies`

### `Deprecated`

## v1.2.5 nf-core/mhcquant "Golden Eagle" - 2019/02/26

### `Added`

### `Fixed`

- sort channels by basename
- fixed psm-level-fdrs

### `Dependencies`

### `Deprecated`

## v1.2.4 nf-core/mhcquant "Golden Eagle" - 2019/02/02

### `Added`

### `Fixed`

- fixed refine_fdr_on_predicted_subset float error

### `Dependencies`

### `Deprecated`

## v1.2.3 nf-core/mhcquant "Golden Eagle" - 2019/02/02

### `Fixed`

- filter out uncommon aminoacids U,X,B,J,Z

## v1.2.2 nf-core/mhcquant "Golden Eagle" - 2019/01/24

### `Added`

### `Fixed`

- default params to false
- change on centroidisation parameter
- small changes on docu

### `Dependencies`

### `Deprecated`

## v1.2.1 nf-core/mhcquant "Golden Eagle" - 2019/01/24

### `Added`

### `Fixed`

- process identical names bug

### `Dependencies`

### `Deprecated`

## v1.2.0 nf-core/mhcquant "Golden Eagle" - 2019/01/19

### `Added`

- Subset FDR refinement option
- Fred2 dependency
- vcf parser and translation to proteins

### `Fixed`

- Documentation

### `Dependencies`

### `Deprecated`

## v1.1.0 nf-core/mhcquant "Black Crow" - 2019/01/04

### `Added`

- optional mhcflurry binding predictions
- peak picking as optional preprocessing step

### `Fixed`

- adapted a few parameters such as the default fdr threshold
- updated documentation

### `Dependencies`

### `Deprecated`

## v1.0.0 nf-core/mhcquant "Naked Chicken" - 2018/11/27

- Initial release of nf-core/mhcquant, created with the [nf-core](http://nf-co.re/) template.
