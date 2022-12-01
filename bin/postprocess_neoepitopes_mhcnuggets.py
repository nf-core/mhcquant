#!/usr/bin/env python
import argparse
import sys
import logging
from collections import defaultdict

# logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("Postprocess Neoepitopes MHCNuggets")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


def main():
    model = argparse.ArgumentParser(description="Postprocess Neoepitopes predicted by MHCNuggets")

    model.add_argument("-i", "--input", type=str, nargs="*", help="predicted class 2 neoepitopes")

    model.add_argument("-n", "--neoepitopes", type=str, help="neoepitopes file")

    args = model.parse_args()

    neoepitope_to_geneID = defaultdict()
    with open(args.neoepitopes, "r") as f:
        # skip header
        f.readline()
        for line in f:
            split = line.split(",")
            neoepitope_to_geneID[split[0]] = split[1]

    for prediction in args.input:
        neoepitope_to_geneID_ic50 = defaultdict()

        with open(prediction, "r") as f:
            # skip header
            f.readline()
            for line in f:
                split = line.split(",")
                neoepitope = split[0].strip()
                geneID = neoepitope_to_geneID[neoepitope].strip()
                ic50 = split[1].strip()
                neoepitope_to_geneID_ic50[neoepitope] = (geneID, ic50)

            with open(prediction + ".csv", "w") as f:
                f.write("peptide,geneID,ic50\n")
                for neoepitope, pair in neoepitope_to_geneID_ic50.items():
                    f.write(neoepitope + "," + pair[0] + "," + pair[1] + "\n")


if __name__ == "__main__":
    main()
