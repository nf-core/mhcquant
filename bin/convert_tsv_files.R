#!/usr/bin/env Rscript
# Written by Marissa Dubbelaar and released under the MIT license.

## Transform the quantified text file to a readable format:
## Arguments:
##   --input         Tabulated data in TSV format which was returned
##                   from the OPENMS_TEXTEXPORTER
##   --outname       Filename for the output repertoire
##   -h  Display help.
## Example: ./convert_tsv_files.R --input input.tsv --outname sample-id.tsv

# Libraries
suppressPackageStartupMessages(library(optparse))
################################################################################
##                       Define commandline arguments                         ##
################################################################################
opt_list <- list(
  make_option(c("--input"), default = NULL,
              help = "Input .tsv file from the openms text exported"),
  make_option(c("--meta_file"), default = NULL,
              help = paste0("Path to the document with meta information ",
                            "about samples")),
  make_option(c("--outname"), default = NULL)
)

opt <- parse_args(OptionParser(option_list=opt_list))
## Read the input content
input <- readLines(opt$input)
meta <- read.table(opt$meta_file, sep = "\t", header = TRUE)
################################################################################
###                           Read the MAP content                           ###
################################################################################
## Get all of the MAP content from the files
map_lines <- grep("^#?MAP", input)
## Obtain the content that starts with (#)MAP
map_content <- input[map_lines]
## Remove the ";" and "#" characters
map_content <- gsub(";|#", "", map_content)
## Make sure that the MAP content can be converted to a dataframe
sample_info <- read.table(text = map_content, sep = "\t", header = TRUE)
## Substitute the identifier that is necessary to perform the matching
identified_samples <- sub("(_aligned)?.mzML", "", sample_info$filename)
## Identify the matching samples based in identifier
matching_samples <-
  as.numeric(sapply(strsplit(identified_samples, "_"), function(x) x[1]))
## Obtain the raw_filenames that have been associated with this sample
raw_filenames <- meta$ReplicateFileName[which(meta$ID == matching_samples)]

## Determine if "msms" in in the raw file names
if (length(grep("msms", raw_filenames)) != 0) {
  ## Add the raw file names to the sample_information
  sample_info$ms_replicate <- sub(".raw", "",
                                  basename(gsub("_", "/", raw_filenames)))
} else {
  ## Otherwise use the ID that is given in the metadata sheet
  sample_info$ms_replicate <- matching_samples
}

################################################################################
###                         Read peptide information                         ###
################################################################################
# Read peptide information
peptide_lines <- grep("^#?PEPTIDE", input)
consensus_lines <- grep("CONSENSUS", input)
if (length(peptide_lines) > 0 && length(consensus_lines) > 0) {
  ## Get all of the peptide information and the belonging consensus information
  peptide_content <-
    paste(input[peptide_lines], input[consensus_lines], sep = "\t")
  peptide_content <- gsub(";|#", "", peptide_content)
  ## Get all of the peptide information and the belonging consensus information
  peptide_info <- read.table(text = peptide_content, sep = "\t", header = TRUE)
  ## Loop through the unique IDs
  for (item in sample_info$id[order(sample_info$id, decreasing = TRUE)]) {
    ## Obtain the filename of the run it belongs to
    filename <- sample_info$ms_replicate[which(sample_info$id == item)]
    ## Substitute the information of the id with the filenames
    colnames(peptide_info) <- gsub(paste0("_", item),
                                   paste0("_", filename),
                                   colnames(peptide_info))
  }
}

################################################################################
###                          Rewrite the MS.* columns                        ###
################################################################################
## Change the "unknown" MS columns for complementing column names
pep_col <- colnames(peptide_info)
## Names are changes using the following REF:
## https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo
pep_col <- gsub("MS.1001491", "Percolator:Q value", pep_col)
pep_col <- gsub("MS.1001492", "Percolator:score", pep_col)
pep_col <- gsub("MS.1001493", "Percolator:PEP", pep_col)
pep_col <- gsub("MS.1002252", "Comet:xcorr", pep_col)
pep_col <- gsub("MS.1002253", "Comet:deltacn", pep_col)
pep_col <- gsub("MS.1002254", "Comet:deltacnstar", pep_col)
pep_col <- gsub("MS.1002255", "Comet:spscore", pep_col)
pep_col <- gsub("MS.1002256", "Comet:sprank", pep_col)
pep_col <- gsub("MS.1002257", "Comet:expectation value", pep_col)
pep_col <- gsub("MS.1002258", "Comet:matched ions", pep_col)
pep_col <- gsub("MS.1002259", "Comet:total ions", pep_col)
## Assign the pep_col variable to the colnames
colnames(peptide_info) <- pep_col

################################################################################
###                         Processing modifications                         ###
################################################################################
## Check the sequences for a modification and trim the ends
patt <- "\\([A-z]*\\)"
mods <- which(grepl(patt, peptide_info$sequence))
current_seq_col <- which(colnames(peptide_info) == "sequence")

peptide_data <- cbind(peptide_info[1:current_seq_col - 1],
                      data.frame(sequence = peptide_info$sequence,
                                 sequence_no_mods =
                                   gsub(patt, "", peptide_info$sequence),
                                 modifications = ""),
                      peptide_info[(current_seq_col + 1):ncol(peptide_info)])

## Loop through the sequences that have modifications and return the content in
## the dataframe
for (hit in mods) {
  ## Obtain the modifications that are present in the different peptides
  out <- gregexpr(patt, peptide_data$sequence[hit])
  ## Obtain the attibutes from the elements in the list with mod matches
  starts <- `attributes<-`(out[[1]], NULL)
  lens <- attr(out[[1]], "match.length")
  ## Output the content as a string
  peptide_data$modifications[hit] <-
    paste(lapply(seq(1, length(starts)), function(n) {
      paste0(substr(peptide_data$sequence[hit], starts[n] - 1, starts[n] - 1),
             starts[n],
             substr(peptide_data$sequence[hit], starts[n],
                    starts[n] + lens[n] - 1))
    }), collapse = ", ")
}

## Save peptide content
write.table(peptide_info, opt$outname, quote = FALSE, sep = "\t",
            row.names = FALSE)
