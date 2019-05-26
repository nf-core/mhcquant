#!/usr/bin/env python
"""
Commandline tool for extracting unique neoepitopes generated from MHCFlurry predictions
and the FRED2 vcf_neoepitope_predictor.py script.

usage: found_neoepitopes.py [-h]
                                -m MZTAB
                                -n NEOEPITOPES
                                -p PREDICTED_PEPTIDES
                               [-f FILEFORMAT {raw, csv, json, pep}]
                                -o OUTPUT

Neoepitope prediction for TargetInsepctor.

optional arguments:
  -h, --help            show this help message and exit
  -m, --mztab           MZTAB
                        Path to the mzab file
  -n, --neoepitopes NEOEPITOPES
                        Path to the neoepitopes input file
  -p, --predicted_peptides PEPTIDES
                        Path to numerous files containing predicted peptides by MHCFlurry
  -f, --file_format {raw, csv, json, pep}
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
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
LOG = logging.getLogger("Found Neoepitopes")
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
    seqs = [l.split()[1] for l in mztab_read if
            l.startswith("PEP") and not 'U' in l and not 'X' in l and not 'Z' in l and not 'J' in l and not 'B' in l]
    seqs_new = list(set(seqs))

    return seqs_new


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
            reader = csv.DictReader(tsvfile, dialect='excel-tab')
            for row in reader:
                neoepitopes.append(row['Sequence'])

        return neoepitopes

    with open(neoepitope_file) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        reader, reader_iter_cp = tee(reader)

        # determine whether any alleles were used for the vcf_neoepitope_predictor script
        # which are not present in the alleles specified for the MHCFlurry prediction
        content = reader_iter_cp.next()
        for allele in alleles:
            try:
                content[allele]
            except KeyError:
                LOG.warning("Allele " + str(allele)
                            + " was specified, but not found in the possible neoepitopes predicted from the vcf file.")
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


def parse_predicted_peptides(predicted_peptides):
    """
    parses the output of MHCFlurry

    :param predicted_peptides: output file of MHCFlurry
    :return: dictionary of alleles to peptide sequences
    """
    HLA_allele_to_peptides = defaultdict(list)

    for index, file in enumerate(predicted_peptides):
        with open(file[0], 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=',')
            current_allele = predicted_peptides[index][1]
            for row in reader:
                HLA_allele_to_peptides[current_allele].append(row['peptide'])

    return HLA_allele_to_peptides


def write_found_neoepitopes(filepath, found_neoepitopes, file_format="csv"):
    """
    writes all unique neoepitopes to a specified file

    :param filepath: path to write the file to
    :param found_neoepitopes: dictionary of alleles to peptide sequences
    :param file_format: json or csv or raw
    """
    # if only list (no bindings)
    if file_format == 'pep':
        with open(filepath + "." + file_format, 'w') as f:
            f.write('\n'.join(str(line) for line in found_neoepitopes))
    elif file_format == "json":
        json.dump(found_neoepitopes + "." + file_format, open(filepath, 'w'))
    elif file_format == "csv":
        dict_list = [found_neoepitopes]
        with open(filepath + "." + file_format, 'w') as f:
            writer = csv.DictWriter(f, dict_list[0].keys())
            writer.writeheader()
            for d in dict_list:
                writer.writerow(d)
    elif file_format == "raw":
        f = open(filepath + "." + file_format, "w")
        f.write(str(found_neoepitopes))
        f.close()
    else:
        LOG.error("Could not write found neoepitopes. Please specify one of the file formats json, csv or raw.")


def main():
    model = argparse.ArgumentParser(description='Neoepitope prediction for TargetInspector.')

    model.add_argument(
        '-p', '--predicted_peptides',
        nargs='+',
        type=str,
        help='One or more allele files containing predicted protein sequences'
    )

    model.add_argument(
        '-n', '--neoepitopes',
        type=str,
        help='All possible predicted neoepitopes'
    )

    model.add_argument(
        '-f', '--file_format',
        type=str,
        default="csv",
        help='File format for output file'
    )

    model.add_argument(
        '-m', '--mztab',
        type=str,
        help='Path to mztab file'
    )

    model.add_argument(
        '-bind', '--predict_bindings',
        action="store_true",
        help='Whether bindings were predicted'
    )

    model.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output file path'
    )

    args = model.parse_args()

    if args.predict_bindings:
        # input files usually look like this: /path/to/A_01_01_predicted_peptides.csv
        # remove any preceding and subsequent letters from the alleles
        # result looks like: A*01:01
        regex = re.compile("[a-zA-Z][_|*][0-9]{2}[_|:][0-9]{2}")
        alleles = map(lambda x: x[:5].replace('_', ':') + x[5:],
                      map(lambda x: x[:2].replace('_', '*') + x[2:],
                          [regex.search(predicted_peptides_file).group() for predicted_peptides_file in
                           args.predicted_peptides]))

        # get alleles to peptide sequences for all found epitopes and possible neoepitopes
        predicted_vcf_neoepitopes = parse_vcf_neoepitopes(args.neoepitopes, alleles)
        found_epitopes = parse_predicted_peptides(zip(args.predicted_peptides, alleles))

        # build the intersection of all found epitopes and possible neoepitopes
        found_neoepitopes = defaultdict(list)
        for allele, sequences in found_epitopes.items():
            found_neoepitopes[allele] = list(set(sequences) & set(predicted_vcf_neoepitopes[allele]))

        write_found_neoepitopes(args.output, found_neoepitopes, args.file_format)

    else:
        # parse all identified peptides and possible neoepitopes
        predicted_vcf_neoepitopes = parse_vcf_neoepitopes(args.neoepitopes, [])
        identified_peptides = parse_mztab(args.mztab)

        # build the intersection of all found epitopes and possible neoepitopes
        found_neoepitopes = list(set(predicted_vcf_neoepitopes) & set(identified_peptides))

        write_found_neoepitopes(args.output, found_neoepitopes, 'pep')


if __name__ == "__main__":
    sys.exit(main())