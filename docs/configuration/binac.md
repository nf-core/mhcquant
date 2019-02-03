# BINAC

You may use the pipeline with the `-profile binac` switch when starting the pipeline. A typical call could work like this for example
```
nextflow run nf-core/mhcquant --mzmls '*.mzML' --fasta 'SWISSPROT_12_2018.fasta' --alleles 'alleles.tsv' --vcf 'variants.vcf' --include_proteins_from_vcf --run_prediction -profile binac
```

