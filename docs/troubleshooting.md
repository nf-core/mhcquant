# nf-core/openmspeptidequant: Troubleshooting

## Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/openmspeptidequant) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
