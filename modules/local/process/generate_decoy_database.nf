// Import generic module functions
include {  saveFiles } from './functions'

params.options = [:]

process GENERATE_DECOY_DB {
    input:
    tuple val(id), val(Sample), file(fastafile)

    output:
    tuple val("$id"), val("$Sample"), file("${fastafile.baseName}_decoy.fasta")
    
    when:
    !params.skip_decoy_generation
 
    script:

     """
     DecoyDatabase  -in ${fastafile} \\
                    -out ${fastafile.baseName}_decoy.fasta \\
                    -decoy_string DECOY_ \\
                    -decoy_string_position prefix
     """
}