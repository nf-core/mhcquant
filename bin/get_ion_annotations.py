#!/usr/bin/env python
__author__ = 'Jonas Scheid'

from typing import Tuple
from pyopenms import *
import pandas as pd
import re
import numpy as np
import argparse
import logging

def parse_arguments() -> Tuple[argparse.ArgumentParser,argparse.Namespace]:
    """
    Purpose: Parse command line arguments
    Output: parser
            args: Contains all parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="""Ion annotator \n Returns two tsv files (all_ions, matching_ions)
    that retrieve all ions with their respective m/z and intensities as well as the ion annotation of ions matching
    to the peptides theoretical spectrum""")

    parser.add_argument('-i', '--input', required=True, nargs='+',
                        help='List of replicate mzML files of a sample')
    parser.add_argument('-idxml', '--filtered_idXML', required=True, type=str,
                        help='FDR filtered idXML file that contains all peptides of the final results')
    parser.add_argument('-p', '--prefix', required=True, type=str,
                        help='Prefix that will be added to the two output files')
    parser.add_argument('-pc', '--precursor_charge', type=str, default='2:3',
                        help='Precursor charge range')
    parser.add_argument('-fmt', '--fragment_mass_tolerance', type=str, default='0.02',
                        help='The fragment mass tolerance of theoretical and experimental spectrum matching')
    parser.add_argument('-a_ions', '--use_a_ions', action='store_true',
                        help='Add a-ions to the theoretical spectrum generation')
    parser.add_argument('-c_ions', '--use_c_ions', action='store_true',
                        help='Add c-ions to the theoretical spectrum generation')
    parser.add_argument('-x_ions', '--use_x_ions', action='store_true',
                        help='Add x-ions to the theoretical spectrum generation')
    parser.add_argument('-z_ions', '--use_z_ions', action='store_true',
                        help='Add z-ions to the theoretical spectrum generation')
    parser.add_argument('-rpp', '--remove_precursor_peak', action='store_true',
                        help='Do not consider precursor masses in the theoretical spectrum generation')

    args = parser.parse_args()

    return parser, args


def generate_theoretical_spectrum(peptide: PeptideIdentification, args: argparse.Namespace) -> MSSpectrum:
    """
    Purpose: Generate theoretical spectrum based on PeptideIdentification
    Output: Theoretical spectrum of input peptide
    """
    sequence = peptide.getHits()[0].getSequence()
    min_charge = int(args.precursor_charge.split(':')[0])
    # If precursor charge ranges between 2:5, the fragment charges range from 1:4.
    # Get the precursor charge information from the idXML file
    min_fragment_charge, max_fragment_charge = min_charge - 1, max(peptide.getHits()[0].getCharge() - 1, 1)

    # Define parameters for theoretical spectrum generation
    tsg = TheoreticalSpectrumGenerator()
    theo_spectrum = MSSpectrum()
    p = tsg.getParameters()
    p.setValue("add_y_ions", "true")
    p.setValue("add_b_ions", "true")
    if args.use_a_ions:
        p.setValue("add_a_ions", "true")
    if args.use_c_ions:
        p.setValue("add_c_ions", "true")
    if args.use_x_ions:
        p.setValue("add_x_ions", "true")
    if args.use_z_ions:
        p.setValue("add_z_ions", "true")
    p.setValue("add_metainfo", "true")
    p.setValue("add_losses", "true")
    p.setValue("add_precursor_peaks", "true")
    p.setValue("add_first_prefix_ion", "true")
    p.setValue("add_abundant_immonium_ions", "true")
    if not args.remove_precursor_peak:
        p.setValue("add_all_precursor_charges", "true")

    # Generate theoretical spectrum
    tsg.setParameters(p)
    tsg.getSpectrum(theo_spectrum, sequence, min_fragment_charge, max_fragment_charge)

    return theo_spectrum


def flatten(ls: list) -> list:
    """
    Purpose: Flatten list of lists into a flat lat
    Output: List
    """
    return [item for sublist in ls for item in sublist]


def __main__():
    parser, args = parse_arguments()
    protein_ids = []
    peptide_ids = []
    # Each peptide ID should only contain one hit in the FDR_filtered idXML
    IdXMLFile().load(args.filtered_idXML, protein_ids, peptide_ids)
    # Define empty lists that collect all the necessary information, which is comprised in DataFrames
    ions = []
    spectra_mz = []
    spectra_intensities = []
    spectra_peptides = []
    is_matching_ion = []
    for i, file in enumerate(args.input):
        # Load the mzML files into the pyopenms structure
        exp = MSExperiment()
        MzMLFile().load(file, exp)
        # Create map containing the spectrum nativeID and the respective scan number
        # TODO: Explore how to access index of spectrum via pyopenms -> DOC is not sufficient to get the necessary information
        native_id_index_map = {spectrum.getNativeID(): index for index, spectrum in enumerate(exp.getSpectra())}
        """native_id_index_map = {spectrum.getNativeID(): int(re.search('scan=(.+)', spectrum.getNativeID()).group(1)) - 1
                               for spectrum in exp.getSpectra()}"""
        # Iterate over the FDR filtered peptideIDs
        for peptide_id in peptide_ids:
            # Access raw spectrum via the spectrum index
            if peptide_id.getMetaValue("spectrum_reference") not in native_id_index_map.keys():
                continue
            spectrum = exp.getSpectrum(native_id_index_map[peptide_id.getMetaValue("spectrum_reference")])
            logging.warn(spectrum.getNativeID())
            logging.warn(peptide_id.getMetaValue("spectrum_reference"))
            # Save mz and intensities of all peaks to comprise them later in a DataFrame
            spectrum_mz, spectrum_intensities = spectrum.get_peaks()
            spectra_mz.append(spectrum_mz)
            spectra_intensities.append(spectrum_intensities)
            sequence = peptide_id.getHits()[0].getSequence()
            spectra_peptides.append(np.repeat(sequence, len(spectrum_mz)))
            is_matching_ion_peptide = np.repeat(False, len(spectrum_mz))
            # Distinguish scans of replicates by checking if the RT of the peptide matches the RT of the spectrum
            # TODO: Come-up with a better solution to distinguish the true scan of the filtered peptide (not via RT)
            if peptide_id.getRT() == round(spectrum.getRT(), 1):
                theo_spectrum = generate_theoretical_spectrum(peptide_id, args)
                alignment = []
                spa = SpectrumAlignment()
                # Specify parameters for the alignment. Pyopenms offers two parameters here
                p = spa.getParameters()
                # Since we look at the MS2 Spectrum we align Da tolerance
                p.setValue("tolerance", int(args.fragment_mass_tolerance))
                p.setValue("is_relative_tolerance", "false")  # false == Da, true == ppm
                spa.setParameters(p)
                # Align both spectra
                spa.getSpectrumAlignment(alignment, theo_spectrum, spectrum)
                # Obtain ion annotations from theoretical spectrum if there are at least two matching ions
                if len(alignment) > 1:
                    for theo_idx, obs_idx in alignment:
                        ion_name = theo_spectrum.getStringDataArrays()[0][theo_idx].decode()
                        ion_charge = theo_spectrum.getIntegerDataArrays()[0][theo_idx]
                        obs_mz = spectrum[obs_idx].getMZ()
                        obs_int = spectrum[obs_idx].getIntensity()
                        ions.append([sequence, ion_name, ion_charge, theo_spectrum[theo_idx].getMZ(), obs_mz, obs_int])

                        is_matching_ion_peptide[obs_idx] = True

            is_matching_ion.append(is_matching_ion_peptide)

    # Save information of matching ions
    matching_ions = pd.DataFrame.from_records(ions, columns=['Peptide', 'Ion_name', 'Ion_charge', 'Theoretical_mass',
                                                             'Experimental_mass', 'Intensity'])
    # Save information of all peaks (including non-matching peaks)
    all_peaks_df = pd.DataFrame([flatten(spectra_peptides), flatten(spectra_mz), flatten(spectra_intensities),
                                 flatten(is_matching_ion)],
                                index=['Peptide', 'Experimental_mass', 'Intensity', 'Is_matching_ion']).transpose()
    all_peaks_df.index.name = 'Ion_ID'
    # Numpy magically converts booleans to floats of 1.0 and 0.0, revert that
    all_peaks_df['Is_matching_ion'] = all_peaks_df['Is_matching_ion'].astype(bool)
    # Add ion ID column to matching ions table for easy reference
    matching_ions['Ion_ID'] = all_peaks_df.index[all_peaks_df['Is_matching_ion']]
    # Write matching ions table
    matching_ions.to_csv(f"{args.prefix}_matching_ions.tsv", sep='\t', index=False)
    all_peaks_df.to_csv(f"{args.prefix}_all_peaks.tsv", sep='\t')


if __name__ == "__main__":
    __main__()
