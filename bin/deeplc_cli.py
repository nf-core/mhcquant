#!/usr/bin/env python
# Written by Jonas Scheid and Steffen Lemke

import click
import logging
import math
import os
import pandas as pd
import sys
import tensorflow as tf
from deeplc import DeepLC
from pyopenms import IdXMLFile
from sklearn.preprocessing import MinMaxScaler

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # Set TensorFlow logging level to suppress warnings
tf.get_logger().setLevel(logging.ERROR)  # Filter out specific warnings

# initate logger
console = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
LOG = logging.getLogger("DeepLC prediction")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)


def parse_idxml(path: str) -> tuple[list, list]:
    """
    Parse idXML file and return PeptideIdentification and ProteinIdentification objects.

    :param path: path to idXML file
    :type path: str
    :return: ProteinIdentification and PeptideIdentification objects
    :rtype: (list, list)
    """
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(path, protein_ids, peptide_ids)

    return protein_ids, peptide_ids


def generate_deeplc_input(peptide_ids: list) -> pd.DataFrame:
    """
    Generate input for DeepLC from PeptideIdentification objects.

    :param peptide_ids: list of PeptideIdentification objects
    :type peptide_ids: list
    :return: Pandas DataFrame containing the input for DeepLC
    :rtype: pd.DataFrame

    """
    data = []
    for peptide_id in peptide_ids:
        tr = peptide_id.getRT()
        for hit in peptide_id.getHits():
            sequence = hit.getSequence()
            unmodified_sequence = sequence.toUnmodifiedString()
            x_corr = hit.getMetaValue("MS:1002252")
            target_decoy = hit.getMetaValue("target_decoy")

            # get all modificaitons
            hit_mods = []
            for pos in range(0, sequence.size()):
                residue = sequence.getResidue(pos)
                if residue.isModified():
                    hit_mods.append("|".join([str(pos + 1), residue.getModificationName()]))
            if hit_mods == []:
                modifications = ""
            else:
                modifications = "|".join(hit_mods)

            data.append([unmodified_sequence, modifications, tr, x_corr, target_decoy])

    df_deeplc_input = pd.DataFrame(data, columns=["seq", "modifications", "tr", "x_corr", "target_decoy"])

    return df_deeplc_input


def generate_calibration_df(df: pd.DataFrame, num_bins: int) -> pd.DataFrame:
    """
    Generates a pandas DataFrame containing calibration peptides for DeepLC.
    The input DataFrame is sorted by measured retention time and sliced into
    bins of equal peptide count. For each bin the peptides with the highest
    x_correlation is selected and returned in a Pandas DataFrame

    :param df: Input DataFrame with retention time of each peptide and xcorr score
    :type df: pd.DataFrame
    :param num_bins: Number of bins/number of calibratoin peptides to be extracted
    :type num_bins: int
    :return: Pandas DataFrame containing calibration peptides equal to index-based num_bins
    :rtype: pd.DataFrame

    """
    # remove decoys
    df = df[df['target_decoy'] != 'decoy']

    # Compute the bin size based on the number of bins
    bin_size = len(df) // num_bins

    # Sort the dataframe by tr values
    sorted_df = df.sort_values('tr')

    # Rows for dataframe
    filtered_row = []

    # Iterate over the bins
    for i in range(num_bins):
        # Get the start and end indices of the current bin
        start_index = i * bin_size
        end_index = start_index + bin_size

        # Get the subset of the dataframe for the current bin
        bin_df = sorted_df.iloc[start_index:end_index]

        # Find the row with the maximum x_corr value in the current bin
        max_row = bin_df.loc[bin_df['x_corr'].idxmax()]

        # Append the max row to the filtered dataframe
        filtered_row.append(max_row)

    # Create DataFrame
    calibration_df = pd.DataFrame(filtered_row)

    return calibration_df.copy()


def generate_calibration_df_with_RT_bins(df: pd.DataFrame, num_bins: int) -> pd.DataFrame:
    """
    Generates a pandas DataFrame containing calibration peptides for DeepLC.
    The input DataFrame is sorted by measured retention time and sliced into bins of equal retention time.
    For each bin the peptides with the highest x_correlation is selected and return in a Pandas DataFrame

    :param df: Input DataFrame with retention time of each peptide and xcorr score
    :type df: pd.DataFrame
    :param num_bins: Number of bins/number of calibratoin peptides to be extracted
    :type num_bins: int
    :return: Pandas DataFrame containing calibration peptides equal to RT-based num_bins
    :rtype: pd.DataFrame
    """
    # remove decoys
    df = df[df['target_decoy'] != 'decoy']

    # Sort the dataframe by tr values
    sorted_df = df.sort_values('tr')

    # Create list of linear bins between min and max tr with num_bins and access dataframe with index
    bin_size = (sorted_df['tr'].max() - sorted_df['tr'].min()) / num_bins

    # Rows for dataframe
    filtered_row = []

    # Iterate over the bins
    for i in range(num_bins):
        # Get the start and end indices of the current bin
        start_tr = sorted_df['tr'].min() + i * bin_size
        end_tr = start_tr + bin_size

        # Get the subset of the dataframe for the current bin
        bin_df = sorted_df[(sorted_df['tr'] >= start_tr) & (sorted_df['tr'] < end_tr)]

        # Find the row with the maximum x_corr value in the current bin
        max_row = bin_df.loc[bin_df['x_corr'].idxmax()]

        # Append the max row to the filtered dataframe
        filtered_row.append(max_row)

    # Create DataFrame
    calibration_df = pd.DataFrame(filtered_row)

    return calibration_df


def min_max_scaler(df: pd.DataFrame) -> pd.DataFrame:
    """
    Scales the predicted retention time values of the input DataFrame to the range of the measured retention time values

    :param df: Input DataFrame with predicted retention time values
    :type df: pd.DataFrame
    :return: DataFrame with scaled predicted retention time values
    :rtype: pd.DataFrame
    """
    scaler = MinMaxScaler((min(df["tr"]), max(df["tr"])))
    df['predicted_RT'] = scaler.fit_transform(df[["predicted_RT"]])

    return df


def run_deeplc(df: pd.DataFrame, calibration_df: pd.DataFrame = None) -> pd.DataFrame:
    dlc = DeepLC()
    if calibration_df is not None:
        dlc.calibrate_preds(seq_df=calibration_df)
        preds = dlc.make_preds(seq_df=df)
        df["predicted_RT"] = preds
    else:
        preds = dlc.make_preds(seq_df=df, calibrate=False)
        df["predicted_RT"] = preds
        df = min_max_scaler(df)

    return df


def add_rt_error(peptide_ids: list, prediction_dict: dict, add_abs_rt_error: bool=False, add_sqr_rt_error: bool=False, add_log_rt_error:bool=False) -> list:
    """
    Adds the error of the predicted retention time in comparison to the measured retention time to each peptide hit.
    Different error scores can be selected.

    :param peptide_ids: list of PeptideIdentification objects
    :type peptide_ids: list
    :param prediction_dict: dictionary containing the predicted retention time for each peptide sequence
    :type prediction_dict: dict
    :param add_abs_rt_error: add absolute RT prediction errors to idXML
    :type add_abs_rt_error: bool
    :param add_sqr_rt_error: add squared RT prediction errors to idXML
    :type add_sqr_rt_error: bool
    :param add_log_rt_error: add log RT prediction errors to idXML
    :type add_log_rt_error: bool
    :return: list of PeptideIdentification objects with added error scores
    :rtype: list
    """
    for peptide_id in peptide_ids:
        # Get measured Retention time
        measured_rt = peptide_id.getRT()

        # Initilaize list for edited hits (with added features)
        new_hits = []
        for hit in peptide_id.getHits():
            sequence = hit.getSequence()
            unmodified_sequence = sequence.toUnmodifiedString()

            # Get modifications
            hit_mods = []
            for pos in range(0, sequence.size()):
                residue = sequence.getResidue(pos)
                if residue.isModified():
                    hit_mods.append("|".join([str(pos + 1), residue.getModificationName()]))
            if hit_mods == []:
                modifications = ""
            else:
                modifications = "|".join(hit_mods)

            predicted_rt = prediction_dict[(unmodified_sequence, modifications)]

            # calculate abs error
            if add_abs_rt_error:
                abs_error = abs(measured_rt - predicted_rt)
                hit.setMetaValue("deeplc_abs_error", abs_error)

            # calculate seq error
            if add_sqr_rt_error:
                sqr_error = abs(measured_rt - predicted_rt)**2
                hit.setMetaValue("deeplc_sqr_error", sqr_error)

            # calcultae log error
            if add_log_rt_error:
                log_error = math.log(abs(measured_rt - predicted_rt))
                hit.setMetaValue("deeplc_log_error", log_error)

            new_hits.append(hit)
        peptide_id.setHits(new_hits)

    return peptide_ids


@click.command()
@click.option('-i', '--input', help='input path of idXML', required=True)
@click.option('-o', '--output', help='output path of idXML',
              required=True)
@click.option('--calibration_mode', type=click.Choice(['idx_bin', 'rt_bin', 'min_max']),
              default='rt_bin', help='Calibration method')
@click.option('--calibration_bins', type=click.IntRange(min=2), default=20,
              help='number of bins for calibration')
@click.option('--add_abs_rt_error', is_flag=True,
              help='add absolute RT prediction errors to idXML')
@click.option('--add_sqr_rt_error', is_flag=True,
              help='add squared RT prediction errors to idXML (default if nothing is selected)')
@click.option('--add_log_rt_error', is_flag=True,
              help='add log RT prediction errors to idXML')
@click.option('--debug', is_flag=True, help='Additionally write out calibration file and deeplc output')
def main(input: str,
         output: str,
         calibration_mode: str,
         calibration_bins: int,
         add_abs_rt_error: bool,
         add_sqr_rt_error: bool,
         add_log_rt_error: bool,
         debug: bool):

    # check if at least one error is selected, if not set squared error to true
    num_true = sum([add_abs_rt_error, add_sqr_rt_error, add_log_rt_error])
    if num_true == 0:
        LOG.info("No error calculation was set, falling back to squared error")
        add_sqr_rt_error = True

    LOG.info("Parse idXML")
    protein_ids, peptide_ids = parse_idxml(input)

    if len(peptide_ids) <= calibration_bins:
        LOG.info("Number of peptide hits is smaller than calibration bins. Skipping deeplc prediction.")
        IdXMLFile().store(output, protein_ids, peptide_ids)
        return 0

    LOG.info("Generate DeepLC input")
    df_deeplc_input = generate_deeplc_input(peptide_ids)

    # Run DeepLC
    if calibration_mode == "rt_bin":
        LOG.info("Run DeepLC with RT bin calibration")
        calibration_df = generate_calibration_df_with_RT_bins(df_deeplc_input, calibration_bins)
        if debug:
            calibration_df.to_csv(output + "_calibration.tsv", index=False, sep="\t")
        df_deeplc_output = run_deeplc(df_deeplc_input, calibration_df)
    elif calibration_mode == "idx_bin":
        LOG.info("Run DeepLC with index bin calibration")
        calibration_df = generate_calibration_df(df_deeplc_input, calibration_bins)
        if debug:
            calibration_df.to_csv(output + "_calibration.tsv", index=False, sep="\t")
        df_deeplc_output = run_deeplc(df_deeplc_input, calibration_df)
    elif calibration_mode == "min_max":
        LOG.info("Run DeepLC with min/max calibration")
        df_deeplc_output = run_deeplc(df_deeplc_input)

    if debug:
        df_deeplc_output.to_csv(output + "_deeplc_output.tsv", index=False, sep="\t")

    # Create map containing the predicted retention time for each peptide sequence and modification
    sequence_to_prediction = {}
    for seq, mods, pred_rt in zip(df_deeplc_output['seq'], df_deeplc_output['modifications'], df_deeplc_output['predicted_RT']):
        sequence_to_prediction[(seq, mods)] = pred_rt

    LOG.info("Add error to idXML")
    peptide_ids_pred_RT = add_rt_error(peptide_ids, sequence_to_prediction, add_abs_rt_error, add_sqr_rt_error, add_log_rt_error)

    LOG.info("Write idXML")
    IdXMLFile().store(output, protein_ids, peptide_ids_pred_RT)

    if debug:
        df_deeplc_input.to_csv(output + "_deeplc_input.tsv", index=False, sep="\t")
        if calibration_mode == "rt_bin" or calibration_mode == "idx_bin":
            calibration_df.to_csv(output + "_calibration.tsv", index=False, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main())
