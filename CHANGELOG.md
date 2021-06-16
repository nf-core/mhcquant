# nf-core/mhcquant: Changelog

## v2.0.0 nf-core/mhcquant "Steel Beagle" - 2021/05/20

### `Added`

DSL1 to DSL2 conversion

    - Different processes based on a unique step in the pipeline
    - Inclusion of one sub workflow: refine fdr on predicted subset
    - The process: openms_cometadapter includes commented lines (which could be used as a reference for future module development)
    - MHCquant pipeline is ran from workflows/mhcquant.nf instead of main

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
- remove specifying input as file dirs eg "data/*.mzML"

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
