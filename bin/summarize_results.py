#!/usr/bin/env python

# Written by Julia Graf and released under MIT license.

import pandas as pd
import numpy as np
from pyopenms import AASequence
from argparse import ArgumentParser
import re
from collections import Counter
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# Setup argument parser
parser = ArgumentParser(description='Extract relevant information for the MultiQC report from the TSV file.')
parser.add_argument(
    "--input",
    required=True,
    nargs=1,
    help="The tsv file returned by the OPENMS_TEXTEXPORTER process."
)

parser.add_argument(
    "--out_prefix",
    required=True,
    nargs=1,
    help="The prefix for the output TSV files."
)

parser.add_argument(
    "--quantify",
    action='store_true',
    help="Whether quantification is enabled or not."
)


def comma_separated_list(value):
    return value.split(',')


parser.add_argument(
    "--columns",
    default=None,
    type=comma_separated_list,
    nargs=1,
    help="A comma-separated list of columns to keep."
)


def parse_multiTSV(file_path):
    """Parse the multi-TSV file output of the TEXTEXPORTER process in quantification mode

    Args:
        file_path (str): The file location of the multi-TSV

    Returns:
        DataFrame: Contains all relevant information for further processing
    """
    peptide_rows = []
    consensus_rows = []
    unassigned_peptide_rows = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("PEPTIDE"):
                peptide_rows.append(line.strip().split('\t')[1:])
            elif line.startswith("CONSENSUS"):
                consensus_rows.append(line.strip().split('\t')[1:])
            elif line.startswith("UNASSIGNEDPEPTIDE"):
                unassigned_peptide_rows.append(line.strip().split('\t')[1:])
            # Extract the header information for all TSV files
            elif line.startswith("#PEPTIDE"):
                peptide_cols = line.strip().split('\t')[1:]
            elif line.startswith("#CONSENSUS"):
                consensus_cols = line.strip().split('\t')[1:]
            elif line.startswith("#UNASSIGNEDPEPTIDE"):
                unassigned_peptide_cols = line.strip().split('\t')[1:]

    # Workaround for OpenMS 3.5.0 TextExporter bug (https://github.com/OpenMS/OpenMS/issues/9120):
    # consensusXML export writes a phantom column in data rows between 'end' and 'FFId_category'
    # that is missing from the header. Remove it to realign columns.
    for rows, cols in [(peptide_rows, peptide_cols), (unassigned_peptide_rows, unassigned_peptide_cols)]:
        if rows and len(rows[0]) > len(cols) and 'end' in cols:
            extra_idx = cols.index('end') + 1
            for i, row in enumerate(rows):
                rows[i] = row[:extra_idx] + row[extra_idx + 1:]

    peptide_df = pd.DataFrame(peptide_rows, columns=peptide_cols)
    consensus_df = pd.DataFrame(consensus_rows, columns=consensus_cols)
    # Concatenate CONSENSUS and PEPTIDE columns
    df = pd.concat([peptide_df, consensus_df], axis=1)
    # Add an additional PSM column
    unassigned_peptide_df = pd.DataFrame(unassigned_peptide_rows, columns=unassigned_peptide_cols)
    psm = unassigned_peptide_df.loc[:, "sequence"].value_counts()
    df["psm"] = df["sequence"].map(psm)
    return df

def strip_modifications(seq_with_mods: str) -> str:
    try:
        aa_seq = AASequence.fromString(seq_with_mods)
        return aa_seq.toUnmodifiedString()
    except Exception as e:
        logging.warning(f"Could not parse sequence {seq_with_mods}: {e}")
        return seq_with_mods

def multi_tsv_not_empty(file_path):
    """A populated multi-TSV (quant) export has PEPTIDE/CONSENSUS data rows (no leading #)."""
    with open(file_path) as f:
        for line in f:
            if line.startswith('PEPTIDE\t') or line.startswith('CONSENSUS\t'):
                return True
    return False


def process_file(file, prefix, quantify, keep_cols):
    """Extract all the relevant information and write it to a TSV file for the MultiQC report

       Args:
           file (str): The file location of the TSV file
           prefix (str): The prefix for all the output files
           quantify (bool): Whether quantification mode is enabled
           keep_cols (list): Set of columns to keep from the original set of columns
       """
    if quantify and multi_tsv_not_empty(file):
        data = parse_multiTSV(file)
        n_psms = np.sum(data["psm"])
    else:
        data = pd.read_csv(file, sep='\t')
        # Remove the special character '#' from the first column name
        data.rename(columns={data.columns[0]: data.columns[0].replace('#', '')}, inplace=True)
        n_psms = 0

    # Check if all required columns are present in the DataFrame
    required_columns = ['sequence', 'accessions', 'mz', 'rt', 'score', 'COMET:xcorr']
    missing_columns = set(required_columns) - set(data.columns)
    if data.shape[0] > 0:  # If the DataFrame is not empty
        if missing_columns:
            raise ValueError(f"The following required columns are missing: {missing_columns}")
    else:
        with open(f"{prefix}_general_stats.csv", "w") as f:
            f.write(f"Sample,# Peptides,# Modified Peptides,# Proteins,# PSMs\n")
            f.write(f"{prefix},0,0,0,0\n")
        data.to_csv(f"{prefix}.tsv", sep='\t', index=False)
        return

    # Remove modification information from the sequence column
    data["peptidoform"] = data["sequence"]
    data["sequence"] = data["sequence"].apply(strip_modifications)

    # ---------------------------------
    # Length distribution plot
    # ---------------------------------

    seq_length = data["sequence"].apply(lambda seq: len(seq))
    seq_length = dict(Counter(seq_length))
    with open(f"{prefix}_peptide_length.csv", "w") as f:
        for length, count in seq_length.items():
            f.write(f"{length},{count / len(data.index)}\n")

    # ---------------------------------
    # General statistics
    # ---------------------------------

    n_peptides = len(set(data["sequence"]))
    n_modified_peptides = sum(1 for s in set(data["peptidoform"]) if '(' in s)
    # Split the accession codes and count each protein accession individually
    n_proteins = len(set([protein for entry in data["accessions"] for protein in entry.split(';')]))
    with open(f"{prefix}_general_stats.csv", "w") as f:
        f.write(f"Sample,# Peptides,# Modified Peptides,# Proteins,# PSMs\n")
        f.write(f"{prefix},{n_peptides},{n_modified_peptides},{n_proteins},{n_psms}\n")

    # ---------------------------------
    # Histograms
    # ---------------------------------

    # Rename IM column to ion_mobility for readability
    if "IM" in data.columns:
        data.rename(columns={"IM": "ion_mobility"}, inplace=True)

    histograms = [[data["mz"].astype(float), f"{prefix}_histogram_mz.csv"],
                  [data["rt"].astype(float), f"{prefix}_histogram_rt.csv"],
                  [data["score"].astype(float), f"{prefix}_histogram_scores.csv"]]

    for values, title in histograms:
        hist, bin_edges = np.histogram(values, bins='auto')
        with open(title, "w") as f:
            for i in range(len(bin_edges) - 1):
                bin_midpoint = (bin_edges[i] + bin_edges[i + 1]) / 2
                f.write(f'{bin_midpoint},{hist[i]}\n')

    # Generate ion mobility histogram if IM data is available (e.g., timsTOF)
    if "ion_mobility" in data.columns:
        im_values = data["ion_mobility"].astype(float)
        hist, bin_edges = np.histogram(im_values, bins='auto')
        with open(f"{prefix}_histogram_im.csv", "w") as f:
            for i in range(len(bin_edges) - 1):
                bin_midpoint = (bin_edges[i] + bin_edges[i + 1]) / 2
                f.write(f'{bin_midpoint},{hist[i]}\n')

    # ---------------------------------
    # Box plots
    # ---------------------------------
    data["COMET:xcorr"].astype(float).to_csv(
        f"{prefix}_xcorr_scores.csv", index=False, header=False
    )
    if 'intensity_cf' in data.columns:
        np.log2(data["intensity_cf"].astype(float)).to_csv(
            f"{prefix}_peptide_intensity.csv",
            index=False,
            header=False
        )

    # Add a column with unique protein accessions next to accessions
    data.insert(data.columns.get_loc('accessions') + 1, 'unique_accessions',
                data['accessions'].map(lambda x: ';'.join(dict.fromkeys(x.split(';')))))

    # Filter the columns down to a user-defined subset of columns
    if keep_cols:
        missing_columns = set(keep_cols) - set(data.columns)
        if missing_columns:
            logging.warning(f"The following columns do not exist in the DataFrame: {missing_columns}")
            keep_cols = [col for col in keep_cols if col not in missing_columns]
        # Retain columns matching rt_*, mz_*, intensity_*, and charge_*
        regex_patterns = [r'^rt_', r'^mz_', r'^intensity_', r'^charge_']
        for pattern in regex_patterns:
            keep_cols.extend([col for col in data.columns if re.match(pattern, col)])
        # Always include unique_accessions next to accessions
        if 'accessions' in keep_cols and 'unique_accessions' not in keep_cols:
            idx = keep_cols.index('accessions') + 1
            keep_cols.insert(idx, 'unique_accessions')
        # Remove duplicates while retaining order
        keep_cols = list(dict.fromkeys(keep_cols))
        data = data.loc[:, keep_cols]

    # Round all floating point values to 5 decimal places to ensure nf-test checksum stability is guaranteed
    float_cols = data.select_dtypes(include=['float']).columns
    data.loc[:, float_cols] = data.loc[:, float_cols].round(5)

    data.to_csv(f"{prefix}.tsv", sep='\t', index=False)


def main():
    args = parser.parse_args()
    if args.columns:
        cols = args.columns[0]
    else:
        cols = None
    process_file(args.input[0],
                 args.out_prefix[0],
                 args.quantify,
                 cols)


if __name__ == '__main__':
    main()
