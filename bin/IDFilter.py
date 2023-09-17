#!/usr/bin/env python
# Written by Jonas Scheid under the MIT license

from pyopenms import *
import pandas as pd
import os
import argparse

def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    :return: parsed arguments
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Filter idXML by a given whitelist of peptides.")
    parser.add_argument("--input", required=True, type=str, help="Input idXML file.")
    parser.add_argument("--whitelist", required=True, type=str, help="IdXML file, which peptide IDs are used as whitelist filter.")
    parser.add_argument("--output", required=True, type=str, help="Filtered idXML file.")

    return parser.parse_args()


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


def filter_run(protein_ids, peptide_ids, whitelist) -> tuple[list, list]:
    """
    Filter Protein and PeptideIdentifications of one run by a whitelist of PeptideIdentifications.

    :param protein_ids: ProteinIdentification objects
    :type protein_ids: list
    :param peptide_ids: PeptideIdentification objects
    :type peptide_ids: list
    :param whitelist: PeptideIdentification objects to keep in the run
    :type whitelist: list
    """
    filter = IDFilter()
    ids_to_keep = [peptide_id for peptide_id in peptide_ids for hit in peptide_id.getHits() if hit.getSequence().toString() in whitelist]
    filter.keepPeptidesWithMatchingSequences(peptide_ids, ids_to_keep, ignore_mods=False)
    # We only want to have unique peptide sequences
    filter.keepBestPerPeptide(peptide_ids, ignore_mods=False, ignore_charges=False, nr_best_spectrum=1)
    filter.removeEmptyIdentifications(peptide_ids)
    # We only want to have protein accessions that are referenced by the fdr-filtered peptide hits
    filter.removeUnreferencedProteins(protein_ids, peptide_ids)

    return protein_ids, peptide_ids


def main():
    args = parse_args()

    # Read idXML files of runs
    protein_ids, peptide_ids = parse_idxml(args.input)

    # Read file containing peptides to keep
    whitelist_protein_ids, whitelist_peptide_ids = parse_idxml(args.whitelist)
    # Get string representation of peptide sequences in fdr_filtered_peptides
    whitelist_peptides = [hit.getSequence().toString() for id in whitelist_peptide_ids for hit in id.getHits()]

    # Filter runs for peptides only in the fdr_filtered_peptides list
    protein_id_filtered, peptide_ids_filtered = filter_run(protein_ids, peptide_ids, whitelist_peptides)

    # Write filtered run to idXML file
    IdXMLFile().store(args.output, protein_id_filtered, peptide_ids_filtered)


if __name__ == "__main__":
    main()
