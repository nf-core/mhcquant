# nf-core/mhcquant: Troubleshooting

## Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.

## Pipeline was killed

Especially when running the pipeline with limited resources for example on your home computer with many large files you will most likely run into memory issues. In fact, if one of the tools along the pipeline requires more memory than available its process might be killed during execution.

Nextflow offers an option to limit the memory used for execution (--max_memory '8 GB'), yet some of the processes might need more memory than specified or nextflow is unable to regulate its memory need.

As the Comet proteomic search engine is the most memory and computation intensive step, one can limit it by setting the --spectrum_batch_size option. Here you specify the number of spectra processed at the same time (default: 500). Setting this parameter lower will require less memory however slow down the overall execution time.

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/mhcquant) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
