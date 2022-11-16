#!/usr/bin/env python
import argparse
import logging
import sys

# logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("Preprocess Neoepitopes MHCNuggets")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


def parse_neoepitopes(filepath):
    with open(filepath, "r") as f:
        # parse and remove the header
        lines = f.read().split("\n")[1:]
        # filter any leftover empty lines
        lines = filter(lambda x: len(x) > 0, lines)
        lines = [elem.split(",")[0] for elem in lines]

        return lines


def write_neoepitopes(neoepitopes, filepath):
    with open(filepath, "w") as f:
        for neoepitope in neoepitopes:
            if neoepitope == neoepitopes[-1]:
                f.write(neoepitope)
            else:
                f.write(neoepitope + "\n")


def main():
    model = argparse.ArgumentParser(description="Neoepitope preprocessing for mhcnuggets")

    model.add_argument("-n", "--neoepitopes", type=str, help="neoepitopes input file")

    model.add_argument(
        "-o", "--output", type=str, help="preprocess neoepitope file for subsequent mhcnuggets prediction"
    )

    args = model.parse_args()

    neoepitopes = parse_neoepitopes(args.neoepitopes)
    write_neoepitopes(neoepitopes, args.output)


if __name__ == "__main__":
    main()
