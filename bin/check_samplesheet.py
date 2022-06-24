#!/usr/bin/env python

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/mhcquant samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    ID  Sample  Condition   ReplicateFileName
    1   WT  A   WT_A.raw
    2   WT  B   WT_B.raw
    3   KO  A   KO_A.raw
    4   KO  B   KO_B.raw
    """

    sample_run_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 4
        HEADER = ["ID", "Sample", "Condition", "ReplicateFileName"]
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        print(header)
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format("\t".join(header), "\t".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split("\t")]
            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            ident, sample, condition, filename = lspl[: len(HEADER)]
            # sample, condition, filename = lspl[1: len(HEADER)]

            ## Check replicate entry is integer
            if not ident.isdigit():
                ident = int(ident)

            # identifier = sample + "_" + condition + "_" + ident
            identifier = ident

            for strCon in [sample, condition]:
                if strCon:
                    if strCon.find(" ") != -1:
                        print_error("Group entry contains spaces!", "Line", line)
                else:
                    print_error("Group entry has not been specified!", "Line", line)

            ## Check MS file extension
            if filename:
                if filename.find(" ") != -1:
                    print_error("FastQ file contains spaces!", "Line", line)

                if not filename.endswith(".raw") and not filename.endswith(".mzML"):
                    print_error(
                        "FastQ file does not have extension '.raw' or '.mzML'!",
                        "Line",
                        line,
                    )
                elif filename.endswith(".raw"):
                    sample_info = [sample, condition, filename, "raw"]
                elif filename.lower().endswith(".mzml"):
                    sample_info = [sample, condition, filename, "mzml"]

                ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, strandedness ]]}
                if sample not in sample_run_dict:
                    sample_run_dict[identifier] = [sample_info]
                else:
                    if sample_info in sample_run_dict[identifier]:
                        print_error("Samplesheet contains duplicate rows!", "Line", line)
                    else:
                        sample_run_dict[identifier].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write("\t".join(["ID", "Sample", "Condition", "Filename", "FileExt"]) + "\n")

            for sample in sorted(sample_run_dict.keys()):
                for idx, val in enumerate(sample_run_dict[sample]):
                    fout.write("\t".join([sample] + val) + "\n")

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
