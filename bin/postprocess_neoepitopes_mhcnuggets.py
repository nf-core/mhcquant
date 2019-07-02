#!/usr/bin/env python
import argparse
from collections import defaultdict


def main():
    model = argparse.ArgumentParser(description='Neoepitope prediction for TargetInspector.')

    model.add_argument(
        '-i', '--input',
        type=str,
        nargs='*',
        help='predicted class 2 neoepitopes'
    )

    model.add_argument(
        '-n', '--neoepitopes',
        type=str,
        help='neoepitopes file'
    )

    model.add_argument(
        '-o', '--output',
        type=str,
        help='postprocessed output files'
    )

    args = model.parse_args()

    neoepitope_to_geneID = defaultdict()
    with open(args.neoepitopes, 'r') as f:
        # skip header
        f.readline()
        for line in f:
            split = line.split(',')
            neoepitope_to_geneID[split[0]] = split[1]

    for prediction in args.input:
        neoepitope_to_geneID_ic50 = defaultdict()

        with open(prediction, 'r') as f:
            # skip header
            f.readline()
            for line in f:
                split = line.split(',')
                neoepitope = split[0]
                geneID = neoepitope_to_geneID[neoepitope]
                ic50 = split[1]
                neoepitope_to_geneID_ic50[neoepitope] = (geneID, ic50)

            with open(prediction, 'w') as f:
                f.write('peptide,geneID,ic50\n')
                for neoepitope, pair in neoepitope_to_geneID_ic50.items():
                    f.write(neoepitope + ',' + pair[0] + ',' + pair[1])


if __name__ == '__main__':
    main()