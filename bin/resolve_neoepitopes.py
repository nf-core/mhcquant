#!/usr/bin/env python
"""
Commandline tool for extracting unique neoepitopes from mztab files and the
and the FRED2 vcf_neoepitope_predictor.py script.
usage: found_neoepitopes.py [-h]
                                -m MZTAB
                                -n NEOEPITOPES
                                [-f FILEFORMAT {raw, csv, json}]
                                -o OUTPUT
Neoepitope prediction for TargetInsepctor.
optional arguments:
    -h, --help            show this help message and exit
    -m, --mztab           MZTAB
                        Path to the mzab file
    -n, --neoepitopes NEOEPITOPES
                        Path to the neoepitopes input file
    -f, --file_format {raw, csv, json}
                        File format to report result in
    -o, --output OUTPUT
                        Path to the output file
"""
import argparse
import csv
import json
import logging
import re
import sys
from itertools import tee
from collections import defaultdict

# logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("Resolve Neoepitopes")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


def parse_mztab(identified_peptides_file):
    """
    parses an mztab file and returns all identified proteins
    :param identified_peptides_file: path to the mztab file
    :return: identified proteins
    """
    mztab = open(identified_peptides_file)
    mztab_read = mztab.readlines()
    mztab.close()

    seq_geneIDs = defaultdict(str)
    for line in mztab_read:
        if line.startswith("PEP"):
            content = line.split("\t")
            seq = content[1]
            geneID = content[2]
            if not "U" in seq and not "X" in seq and not "Z" in seq and not "J" in seq and not "B" in seq:
                seq_geneIDs[seq] = geneID

    return seq_geneIDs


def parse_vcf_neoepitopes(neoepitope_file, alleles):
    """
    parses the output of the VCF neoepitope predictor script
    :param neoepitope_file: output file of VCF neoepitope script
    :param alleles: all alleles which were used for the VCF neoepitope script
    :return: dictionary of alleles to peptide sequences
    """
    HLA_allele_to_peptides = defaultdict(list)

    if not alleles:
        neoepitopes = []
        with open(neoepitope_file) as tsvfile:
            reader = csv.DictReader(tsvfile, dialect="excel-tab")
            for row in reader:
                neoepitopes.append(row["Sequence"])

        return neoepitopes

    with open(neoepitope_file) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect="excel-tab")
        reader, reader_iter_cp = tee(reader)

        # determine whether any alleles were used for the vcf_neoepitope_predictor script
        # which are not present in the alleles specified for the MHCFlurry prediction
        content = reader_iter_cp.next()
        for allele in alleles:
            try:
                content[allele]
            except KeyError:
                LOG.warning(
                    "Allele "
                    + str(allele)
                    + " was specified, but not found in the possible neoepitopes predicted from the vcf file."
                )
                alleles.remove(allele)

        # allele with highest affinity score wins
        for row in reader:
            max_val = 0
            max_allele = ""
            for allele in alleles:
                if row[allele] > max_val:
                    max_allele = allele
                    max_val = row[allele]
            HLA_allele_to_peptides[max_allele].append(row["Sequence"])

    return HLA_allele_to_peptides


def write_found_neoepitopes(filepath, found_neoepitopes, file_format="csv"):
    """
    writes all unique neoepitopes to a specified file
    :param filepath: path to write the file to
    :param found_neoepitopes: dictionary of alleles to peptide sequences
    :param file_format: json or csv or raw
    """
    # if only list (no bindings)
    if file_format == "pep":
        with open(filepath + "." + file_format, "w") as f:
            f.write("Peptide sequence" + "\t" + "geneID")
            f.write("\n".join(str(seq) + "\t" + str(geneID) for seq, geneID in found_neoepitopes.items()))
    elif file_format == "json":
        json.dump(found_neoepitopes, open(filepath + "." + file_format, "w"))
    elif file_format == "csv":
        with open(filepath + "." + file_format, "w") as csv_file:
            writer = csv.writer(csv_file)
            header = ["Peptide sequence", "geneID"]
            writer.writerow(header)
            for seq, geneID in found_neoepitopes.items():
                writer.writerow([seq, geneID])
    elif file_format == "raw":
        f = open(filepath + "." + file_format, "w")
        f.write(str(found_neoepitopes))
        f.close()
    else:
        LOG.error("Could not write found neoepitopes. Please specify one of the file formats json, csv or raw.")


def main():
    model = argparse.ArgumentParser(
        description="Neoepitope resolvement from mztab and possible vcf determined neoepitopes."
    )

    model.add_argument("-n", "--neoepitopes", type=str, help="All possible predicted neoepitopes")

    model.add_argument("-m", "--mztab", type=str, help="Path to mztab file")

    model.add_argument("-f", "--file_format", type=str, default="csv", help="File format for output file")

    model.add_argument("-o", "--output", type=str, required=True, help="Output file path")

    args = model.parse_args()

    # parse all identified peptides and possible neoepitopes
    predicted_vcf_neoepitopes = parse_vcf_neoepitopes(args.neoepitopes, [])
    identified_peptides_to_geneIDs = parse_mztab(args.mztab)

    # build the intersection of all found epitopes and possible neoepitopes
    found_neoepitopes = list(set(predicted_vcf_neoepitopes) & set(identified_peptides_to_geneIDs.keys()))
    LOG.info(str(len(found_neoepitopes)) + ' Neoepitopes were found. Examine "found_neoepitopes.csv" for details.')
    found_neoepitopes_to_geneIDs = {k: v for k, v in identified_peptides_to_geneIDs.items() if k in found_neoepitopes}

    write_found_neoepitopes(args.output, found_neoepitopes_to_geneIDs, args.file_format)


if __name__ == "__main__":
    sys.exit(main())
