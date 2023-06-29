#!/usr/bin/env python
# Written by Jonas Scheid and Steffen Lemke

import click
import logging
import numpy as np
import pandas as pd
import sys
from ms2pip.ms2pipC import MS2PIP
from pyopenms import IdXMLFile, ModificationsDB


# initate logger
console = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("MS2pip prediction")
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


def peptide_ids_to_peprec_dataframe(peptide_ids: list, hit_idx: int = 0) -> pd.DataFrame:
    """
    All the peptide identifications are parsed into a DataFrame in the style of
    a PEPREC file (https://github.com/compomics/ms2pip#peprec-file).

    :param peptide_ids: List containing PeptideIdentification
    :type peptide_ids: list
    :param hit_idx: hit index to generate a peprec
    :type hit_idx: int
    :return: peprec pandas dataframe
    :rtype: pd.DataFrame
    """

    columns = ["spec_id", "modifications", "peptide", "charge"]
    data = []
    spectrum_reference_to_seq = {}

    for peptide_id in peptide_ids:
        if len(peptide_id.getHits()) <= hit_idx:
            continue
        hit = peptide_id.getHits()[hit_idx]
        spectrum_reference = peptide_id.getMetaValue("spectrum_reference")

        charge = hit.getCharge()
        sequence = hit.getSequence()
        unmodified_sequence = sequence.toUnmodifiedString()

        spectrum_reference_to_seq[spectrum_reference] = str(sequence)

        hit_mods = []
        for pos in range(0, sequence.size()):
            residue = sequence.getResidue(pos)
            if residue.isModified():
                hit_mods.append("|".join([str(pos + 1), residue.getModificationName()]))
        if hit_mods == []:
            modifications = "-"
        else:
            modifications = "|".join(hit_mods)

        data.append([spectrum_reference, modifications, unmodified_sequence, charge])

    return pd.DataFrame(data, columns=columns), spectrum_reference_to_seq


def get_complete_spectrum_correlation(df_ms2pip_output: pd.DataFrame, method: str) -> pd.DataFrame:
    """
    Get correlation coefficient for each predicted spectrum vs the measured one

    :param df_ms2pip_output: pandas dataframe of the ms2pip output with individual ion prediction values and DeepLC RT prediction
    :type hit_idx: pd.DataFrame
    :return: dict {<scan_nr>: <Pearson correlation coefficient>, <scan_nr>: {}, ... }
    :rtype: pd.DataFrame
    """
    scannr_to_total_corr = {}
    grouped_spec = df_ms2pip_output.groupby("spec_id")
    correlations_spec = grouped_spec[["prediction", "target"]].corr(method=method)

    for group, corr in correlations_spec.groupby(level=[0, 1]):
        correlation_value = corr.iloc[0, 1]
        spec_id = group[0]
        if group[1] == "prediction":
            if np.isnan(correlation_value):
                correlation_value = 0
            scannr_to_total_corr[spec_id] = correlation_value

    data = {
        "ScanNr": scannr_to_total_corr.keys(),
        "ion_corr": scannr_to_total_corr.values(),
    }
    df = pd.DataFrame.from_dict(data)

    return df


def generate_params_config(
    fixed_modifications: list,
    variable_modifications: list,
    model_name: str,
    fragment_error: float,
) -> dict:
    """
    Generate the MS2PIP configuration file.

    :param fixed_modifications: List of fixed modifications to consider
    :type fixed_modifications: list
    :param modifications: List of modifications to consider
    :type modifications: list
    :param model_name: Name of the model to use
    :type model_name: str
    :param fragment_error: Fragment error to use
    :type fragment_error: float
    :return: MS2PIP configuration file
    :rtype: dict
    """
    mods = set(fixed_modifications.strip().split(",") + variable_modifications.strip().split(","))
    # Remove empty strings
    mods = [mod for mod in mods if mod]
    params = {
        "ms2pip": {
            "ptm": [
                f"{ModificationsDB().getModification(mod).getId()},{ModificationsDB().getModification(mod).getDiffMonoMass()},opt,{ModificationsDB().getModification(mod).getOrigin()}"
                for mod in mods
            ],
            "model": model_name,
            "frag_error": fragment_error,
            "out": "csv",
            "sptm": [],
            "gptm": [],
        }
    }
    return params


@click.command()
@click.option("--input_idxml", help="input path of idXML", required=True)
@click.option("--input_mzml", help="input path of mzML", required=True)
@click.option("--output_idxml", help="output path of idXML", required=True)
@click.option(
    "--num_hits",
    type=click.IntRange(min=1),
    default=1,
    help="number of peptides hits",
)
@click.option(
    "--model_name",
    type=str,
    help="Name of MS2pip model (https://github.com/compomics/ms2pip#specialized-prediction-models)",
)
@click.option(
    "--model_path",
    type=str,
    help="path to MS2pip model",
)
@click.option(
    "--fragment_error",
    type=float,
    help="Fragment mass error in Da",
)
@click.option(
    "--variable_mods",
    type=str,
    help="List of variable modifications",
)
@click.option(
    "--fixed_mods",
    type=str,
    help="List of fixed modifications",
)
@click.option("--add_pearson", is_flag=True, help="add pearson spectrum simliartity")
@click.option(
    "--num_cpus",
    type=int,
    help="number of cpus to use",
)
def main(
    input_idxml: str,
    input_mzml: str,
    output_idxml: str,
    num_hits: int,
    model_name: str,
    model_path: str,
    fragment_error: float,
    variable_mods: str,
    fixed_mods: str,
    add_pearson: bool,
    num_cpus: int,
):
    LOG.info("Parse idXML")
    protein_ids, peptide_ids = parse_idxml(input_idxml)

    LOG.info("Generate params file for MS2pip")
    params = generate_params_config(fixed_mods, variable_mods, model_name, fragment_error)

    LOG.info("Make MS2pip predictions")
    scan_nr_seq_to_corr = {}
    for hit_idx in range(num_hits):  # number of hits to consider
        df_peprec, scan_nr_to_seq = peptide_ids_to_peprec_dataframe(peptide_ids, hit_idx)
        ms2pip = MS2PIP(pep_file=df_peprec, spec_file=input_mzml, params=params, return_results=True, num_cpu=num_cpus)
        predictions = ms2pip.run()
        correlation_df = get_complete_spectrum_correlation(predictions, "pearson")

        for scan_nr, ion_corr in zip(correlation_df["ScanNr"], correlation_df["ion_corr"]):
            sequence = scan_nr_to_seq[scan_nr]
            scan_nr_seq_to_corr[(scan_nr, sequence)] = ion_corr

    LOG.info("Add correlations scores to peptide identifications")
    for peptide_id in peptide_ids:
        spectrum_reference = peptide_id.getMetaValue("spectrum_reference")
        new_hits = []
        for hit in peptide_id.getHits():
            sequence = str(hit.getSequence())
            if (spectrum_reference, sequence) in scan_nr_seq_to_corr.keys():
                hit.setMetaValue(
                    "spectrum_correlation",
                    scan_nr_seq_to_corr[(spectrum_reference, sequence)],
                )
            else:
                LOG.info(f"No correlation could be computed for {str(sequence)}")
                hit.setMetaValue("spectrum_correlation", 0)
            new_hits.append(hit)
        peptide_id.setHits(new_hits)

    LOG.info("Write idXML")
    IdXMLFile().store(output_idxml, protein_ids, peptide_ids)

    return 0


if __name__ == "__main__":
    sys.exit(main())
