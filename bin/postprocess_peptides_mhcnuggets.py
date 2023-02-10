#!/usr/bin/env python
import argparse
import sys
import logging
from collections import defaultdict

# logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("Postprocess Peptides MHCNuggets")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


def main():
    model = argparse.ArgumentParser(description="Postprocess Neoepitopes predicted by MHCNuggets")

    model.add_argument("-i", "--input", type=str, nargs="*", help="predicted class 2 peptides")

    model.add_argument("-p", "--peptides_seq_ID", type=str, help="peptides to seq_ID csv file")

    model.add_argument("-o", "--output", type=str, help="output file name")

    args = model.parse_args()

    peptide_to_geneID = defaultdict()
    with open(args.peptides_seq_ID, "r") as f:
        for line in f:
            split = line.split(",")
            peptide_to_geneID[split[0]] = split[1]

    for prediction in args.input:
        peptide_to_geneID_ic50 = defaultdict()

        with open(prediction, "r") as f:
            # skip header
            f.readline()
            for line in f:
                split = line.split(",")
                peptide = split[0].strip()
                geneID = peptide_to_geneID[peptide].strip()
                ic50 = split[1].strip()
                peptide_to_geneID_ic50[peptide] = (geneID, ic50)

            with open(args.output, "w") as f:
                f.write("peptide,geneID,ic50\n")
                for peptide, pair in peptide_to_geneID_ic50.items():
                    f.write(peptide + "," + pair[0] + "," + pair[1] + "\n")


if __name__ == "__main__":
    main()
