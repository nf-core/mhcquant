# nf-core/mhcquant: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.7.0dev - [date]

### `Added`

- Added `PYOPENMS_CHROMATOGRAMEXTRACTOR` extracting MS1 Chromatograms and visualize them in multiQC report [#329](https://github.com/nf-core/mhcquant/pull/329)
- Added `OPENMS_IDMASSACCURACY` and `DATAMASH_HISTOGRAM` to compute fragment mass errors and visualizte them in multiQC report [#332](https://github.com/nf-core/mhcquant/pull/332)
- Added global fdr evaluation in new local subworkflow `RESCORE` [#338](https://github.com/nf-core/mhcquant/pull/338)
- Added `-weights` parameter in `OPENMS_PERCOLATORADAPTER` and visualize the median feature weights in multiQC report [#347](https://github.com/nf-core/mhcquant/pull/347)
- Added flag `generate_speclib` that will generate a spectrum library for DIA searches with EasyPQP [#349](https://github.com/nf-core/mhcquant/pull/349)
- Replace local with nf-core modules [#350](https://github.com/nf-core/mhcquant/pull/347)

### `Fixed`

- Template update 3.0.2 [#337](https://github.com/nf-core/mhcquant/pull/337)
- Template update 3.1.1 [#346](https://github.com/nf-core/mhcquant/pull/346)
- Template update 3.1.2 [#354](https://github.com/nf-core/mhcquant/pull/354)
- Template update 3.2.0 [#356](https://github.com/nf-core/mhcquant/pull/356)
- Bump OpenMS version 3.1.0 -> 3.3.0 [#358](https://github.com/nf-core/mhcquant/pull/358)

## v2.6.0 - nfcore/mhcquant "Mr Bob" - 2024/06/17

### `Added`

- Added MS²Rescore module with the underlying python CLI [#293](https://github.com/nf-core/mhcquant/pull/293)
- Added support for handling various archive formats: `d|d.tar.gz|d.tar|d.zip|mzML.gz|raw|RAW|mzML` [#323](https://github.com/nf-core/mhcquant/pull/323)
- Added test for timsTOF data [#323](https://github.com/nf-core/mhcquant/pull/323)
- Added new flag `ms2pip_model_dir`, which allows specifying a cache directory for ms2pip models [#322](https://github.com/nf-core/mhcquant/pull/322)

### `Fixed`

- Create only one decoy database [#299](https://github.com/nf-core/mhcquant/pull/299)
- Template update 2.11 [#300](https://github.com/nf-core/mhcquant/pull/300)
- Template update 2.12 [#303](https://github.com/nf-core/mhcquant/pull/303)
- Use `groupKey` to streamline group-wise processing [#310](https://github.com/nf-core/mhcquant/pull/310)
- Replace `PYOPENMS_IDFILTER` with `OPENMS_IDFILTER` [#310](https://github.com/nf-core/mhcquant/pull/310)
- Added nf-core modules [#310](https://github.com/nf-core/mhcquant/pull/310)
- Template update 2.13 [#311](https://github.com/nf-core/mhcquant/pull/311)
- Template update 2.13.1 [#313](https://github.com/nf-core/mhcquant/pull/313)
- Template update 2.14.1 [#320](https://github.com/nf-core/mhcquant/pull/320)
- Added stubs to local modules [#326](https://github.com/nf-core/mhcquant/pull/326)

### `Changed`

- Set identification mode as default and rename `--skip_quantification` to `--quantify` [#323](https://github.com/nf-core/mhcquant/pull/323)

### `Deprecated`

- Removed MS²PIP and DeepLC modules. These feature generators are now called via the MS²Rescore framework [#293](https://github.com/nf-core/mhcquant/pull/293)

## v2.5.0 - nfcore/mhcquant "Angry Bird" - 2023/10/09

### `Added`

- Support for brukers tdf format by adding tdf2mzml converter [#263](https://github.com/nf-core/mhcquant/issues/263)
- DeepLC retention time prediction
- MS2PIP peak intensity prediction
- Added OpenMS FileFilter to clean mzml after parsing to remove artifacts like empty spectra or precursors with charge 0 (optional)
- Made file extension check case insensitive
- Added option to provide a default comet parameters file
- Optimize resource allocations
- Template update 2.9 [#274](https://github.com/nf-core/mhcquant/pull/274)
- Improved quantification such that merged FDR-filtered runs can be quantified properly
- Template update 2.10 [#282](https://github.com/nf-core/mhcquant/pull/282)

### `Fixed`

- [#266](https://github.com/nf-core/mhcquant/pull/266) New OpenMS version 3.0.0 fixes duplicated ID bug [#250](https://github.com/nf-core/mhcquant/issues/250)

### `Dependencies`

- [#266](https://github.com/nf-core/mhcquant/pull/266) Switched from OpenMS version 2.8.0 to newest version 3.0.0 [#265](https://github.com/nf-core/mhcquant/issues/265)
- [#266](https://github.com/nf-core/mhcquant/pull/266) Bumped ThermoRawFileParser version from 1.4.0 to 1.4.2

### `Deprecated`

- OpenMS RT prediction

## v2.4.1 nfcore/mhcquant "Young Shark" (patch) - 2023/04/04

### `Added`

- Added low resolution settings (e.g. Iontrap) [#254](https://github.com/nf-core/mhcquant/pull/254)

### `Fixed`

- Increased comet search, through altering the spectrum_batch_size from 500 to 0
- [#249](https://github.com/nf-core/mhcquant/pull/249) - nf-core template update (version 2.7.2)
- [#258](https://github.com/nf-core/mhcquant/pull/258) - Adjusted decoy strategy to reverse [#255](https://github.com/nf-core/mhcquant/issues/255) and made consistent fdr-level flags [#228](https://github.com/nf-core/mhcquant/issues/228)
- [#845](https://github.com/nf-core/test-datasets/pull/845) - Adjusted nf-core test data set [#233](https://github.com/nf-core/mhcquant/issues/233)

### `Dependencies`

### `Deprecated`

## v2.4.0 nfcore/mhcquant "Maroon Gold Boxer" - 2022/12/02

Initial release of nf-core/mhcquant, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Additional output from `CometAdapter` (generated with the parameter `--pin_out`)
- Folder structure within the `intermediate_results` folder to retrace the outcome files easier
- `OPENMS_FALSEDISCOVERYRATE` and `OPENMS_IDFILTER_FOR_ALIGNMENT` are now included in the first quantification step
- Altered the outcome content with the inclusion of the different folder structure
- Updated the mhcquant_web.png in the `assets` folder
- [#229](https://github.com/nf-core/mhcquant/pull/229) Add ion annotation feature requested in [#220](https://github.com/nf-core/mhcquant/issues/220)
- [#235](https://github.com/nf-core/mhcquant/issues/235) Add the `annotate_ions` parameter to enable/disable the ion annotation feature (default is false)

### `Fixed`

- Resolved issues with `SAMPLESHEET_CHECK`
- Fix for the `peakpickerhires`, mzml files generated from input raw files are now seen as input for this step as well
- `PRE_QUANTIFICATION` is renamed to `MAP_ALIGNMENT` to indicate that the alignment (and the complementing processes) of the different maps happens here
- `POST_QUANTIFICATION` is renamed to `PROCESS_FEATURE` since the feature identification and processing is done here
- Outcome of `OPENMS_FEATUREFINDERIDENTIFICATION` got lost during one of the previous updates, this is reintroduced
- `OPENMS_TEXTEXPORTER_UNQUANTIFIED` and `OPENMS_TEXTEXPORTER_QUANTIFIED` return only significant hits again
- [#226](https://github.com/nf-core/mhcquant/pull/226) - nf-core template update (version 2.6)
- [#230](https://github.com/nf-core/mhcquant/issues/230) - Issue with `OPENMS_MZTABEXPORTER_QUANT`
- [#236](https://github.com/nf-core/mhcquant/issues/236) - Resolved issue with `PYOPENMS_IONANNOTATOR`
- Fix for an inconsistent mzml channel issue
- [#241](https://github.com/nf-core/mhcquant/issues/241) - Fix of the HLA allele annotation in the help of the `allele_sheet` parameter

### `Dependencies`

- Updated the multiQC module

| Dependency            | Old version | New version |
| --------------------- | ----------- | ----------- |
| `MultiQC`             | 1.11        | 1.12        |
| `OpenMS`              | 2.6.0       | 2.8.0       |
| `OpenMS thirdparty`   | 2.6.0       | 2.8.0       |
| `pyOpenMS`            | -           | 2.8         |
| `thermorawfileparser` | 1.3.4       | 1.4.0       |

### `Deprecated`

- `OPENMS_TEXTEXPORTER_PSMS` was removed due to the outcome of the comet adapter step

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
- Adjustments of the `PRE_QUANTIFICATION` subworkflow: `OPENMS_FALSEDISCOVERYRATE`, `OPENMS_IDFILTER_FOR_ALIGNMENT`, and `OPENMS_TEXTEXPORTER_SINGLE`
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
