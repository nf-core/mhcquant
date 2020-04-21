# nf-core/mhcquant: Changelog

## v1.5 nf-core/mhcquat "Flying Fish" - 2020/04/18

### `Fixed`

- raise OpenMS version to 2.5
- adapt workflow accoringly with new options
- integrate sample, allele and vcf sheets instead of file dirs
- remove specifying input as file dirs eg "data/*.mzML"
- branched mzML/raw input
- introduce param to skip quantification

### `Dependencies`

### `Deprecated`

## v1.4 nf-core/mhcquat "Blue Kingfisher" - 2020/03/18

### `Fixed`

- output file order in intermediate results
- increased run times for MS search and variant translation
- Raw File Reading
- RT prediction
- Quantification FDR
- Variant pass filter
- nf-core template update 1.8 and 1.9
- Added version numbers of mhcnuggets and Fred2

### `Dependencies`

### `Deprecated`

## v1.3 nf-core/mhcquant "Red Parrot" - 2019/08/03

### `Fixed`

- nf-core template update
- x,z,a,c ions
- quantification fdr
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

## v1.2.6 nf-core/mhcquant "Golden Eagle" - 2019/03/05

### `Fixed`

- linear retention time alignment
- refine_fdr README

## v1.2.5 nf-core/mhcquant "Golden Eagle" - 2019/02/26

### `Fixed`

- sort channels by basename
- fixed psm-level-fdrs

## v1.2.4 nf-core/mhcquant "Golden Eagle" - 2019/02/02

### `Fixed`

- fixed refine_fdr_on_predicted_subset float error

## v1.2.3 nf-core/mhcquant "Golden Eagle" - 2019/02/02

### `Fixed`

- filter out uncommon aminoacids U,X,B,J,Z

## v1.2.2 nf-core/mhcquant "Golden Eagle" - 2019/01/24

### `Fixed`

- default params to false
- change on centroidisation parameter
- small changes on docu

## v1.2.1 nf-core/mhcquant "Golden Eagle" - 2019/01/24

### `Fixed`

- process identical names bug

## v1.2.0 nf-core/mhcquant "Golden Eagle" - 2019/01/19

### `Added`

- Subset FDR refinement option
- Fred2 dependency
- vcf parser and translation to proteins

### `Fixed`

- Documentation

## v1.1.0 nf-core/mhcquant "Black Crow" - 2019/01/04

### `Added`

- optional mhcflurry binding predictions
- peak picking as optional preprocessing step

### `Fixed`

- adapted a few parameters such as the default fdr threshold
- updated documentation

## v1.0.0 nf-core/mhcquant "Naked Chicken" - 2018/11/27

- Initial release of nf-core/mhcquant, created with the [nf-core](http://nf-co.re/) template.
