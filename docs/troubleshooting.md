# nf-core/mhcquant: Troubleshooting

## Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.

## Pipeline was killed

Especially when running the pipeline with limited resources for example on your home computer with many large files you will most likely run into memory issues. In fact, if one of the tools along the pipeline requires more memory than available its process might be killed during execution.

Nextflow offers an option to limit the memory used for execution (--max_memory '8 GB'), yet some of the processes might need more memory than specified or nextflow is unable to regulate its memory need.

As the Comet proteomic search engine is the most memory and computation intensive step, one can limit it by setting the --spectrum_batch_size option. Here you specify the number of spectra processed at the same time (default: 500). Setting this parameter lower will require less memory however slow down the overall execution time.

## Variant calling file format problems

Variant calling files come in very different formats and annotation styles. While we tried to support most common formats it could be that our integrated vcf_reader can not parse every format correctly resulting in an error. As a general guideline we support VCF format < 4.2 and annotation styles of VEP, ANNOVAR and SNPEFF. Moreover, insertions and deletions might not be parsed properly and one possibility could be to focus on SNPs only in your data. A last work around could be to translate variants into a protein database independantly of MHCquant and use this as input for the database search.

## Fasta database problems

Make sure your fasta database does not contain empty lines or special characters in Sequences.

## Issues with comet.exe on arch based systems

We noticed that comet.exe sometimes crashes on arch based systems using the docker or singularity profiles. It is possible that this is due to the docker package being build with a different compiler than comet and therefore this could potentially lead to header conflicts. Since we are creating the singularity image from the docker image the error is propagated to the singularity profile as well. Possibly changing systems or running MHCquant in a virtual machine is then the only solution.

## Extra resources and getting help

If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/mhcquant) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
