#!/usr/bin/env python
# Written by Jonas Scheid under the MIT license

import logging
import csv
import argparse
import matplotlib.pyplot as plt

import pandas as pd
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
    chromatogram = [(spectrum.getRT() / 60 , spectrum.calculateTIC()) for spectrum in exp.getSpectra() if spectrum.getMSLevel() == 1]
    logging.info(f'Found {len(chromatogram)} MS1 Spectra')
    logging.info(f'RT range: {round(chromatogram[0][0],2)} - {round(chromatogram[-1][0],2)} [min]')
    # Create pandas df
    chromatogram_df = pd.DataFrame(chromatogram, columns=['RT', 'TIC'])
    # bin data into minutes and take the mean of the TIC
    chromatogram_df = chromatogram_df.groupby('RT').mean().reset_index()
    # Add RT=0 and Intensity=0 to start and end of chromatogram_df
    start = pd.DataFrame([{'RT': 0, 'TIC': 0}])
    end = pd.DataFrame([{'RT': chromatogram_df['RT'].max(), 'TIC': 0}])
    # Concatenate the DataFrames
    chromatogram_df = pd.concat([start, chromatogram_df, end], ignore_index=True)

    # Write to csv
    chromatogram_df.to_csv(output_file, index=False, header=False)

if __name__ == '__main__':
    main()
