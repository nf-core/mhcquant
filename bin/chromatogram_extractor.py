#!/usr/bin/env python
# Written by Jonas Scheid under the MIT license

import logging
import csv
import argparse
import matplotlib.pyplot as plt

import pyopenms as oms

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def parse_arguments():
    parser = argparse.ArgumentParser(description='Exctract TICs of MS1 Spectra')
    parser.add_argument('-in','--input', type=str, help='Path to the spectrum file')
    parser.add_argument('-out', '--output', type=str, help='Path to the output CSV file containing RT and TIC of Precursors')
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_file = args.input
    output_file = args.output

    # Load the mzML file
    logging.info(f'Loading file: {input_file}')
    exp = oms.MSExperiment()
    mzml_file = oms.MzMLFile()
    mzml_file.load(input_file, exp)

    # Get RT and Spectrum TIC of MS1 Spectra
    chromatogram = [(spectrum.getRT() / 60, spectrum.calculateTIC()) for spectrum in exp.getSpectra() if spectrum.getMSLevel() == 1]
    logging.info(f'Found {len(chromatogram)} MS1 Spectra')
    logging.info(f'RT range: {round(chromatogram[0][0],2)} - {round(chromatogram[-1][0],2)} [min]')
    # Add (0,0) to start and end of chromatogram
    chromatogram = [(0, 0)] + chromatogram + [(chromatogram[-1][0], 0)]

    # Write to csv
    with open(output_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(chromatogram)

if __name__ == '__main__':
    main()
