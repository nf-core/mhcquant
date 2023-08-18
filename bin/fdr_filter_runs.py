#!/usr/bin/env python
# Written by Jonas Scheid

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
    parser = argparse.ArgumentParser(description="Filter runs by hits in FDR-filtered Percolator output.")
    parser.add_argument("--input", required=True, type=str, help="Input idXML file.")
    parser.add_argument("--pout", required=True, type=str, help="Percolator output file containing fdr filtered hits merged from all replicates.")
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
    filter.keepBestPerPeptide(peptide_ids, ignore_mods=False, ignore_charges=True, nr_best_spectrum=1)
    filter.removeEmptyIdentifications(peptide_ids)
    # We only want to have protein accessions that are referenced by the fdr-filtered peptide hits
    filter.removeUnreferencedProteins(protein_ids, peptide_ids)

    return protein_ids, peptide_ids


def main():
    args = parse_args()

    # Read filtered percolator output
    pout_protein_ids, pout_peptide_ids = parse_idxml(args.pout)

    # Read idXML files of runs
    run_protein_ids, run_peptide_ids = parse_idxml(args.input)

    # Get string representation of peptide sequences in fdr_filtered_peptides
    pout_peptides = [hit.getSequence().toString() for id in pout_peptide_ids for hit in id.getHits()]

    # Filter runs for peptides only in the fdr_filtered_peptides list
    run_protein_id_filtered, run_peptide_ids_filtered = filter_run(run_protein_ids, run_peptide_ids, pout_peptides)

    # Write filtered run to idXML file
    IdXMLFile().store(args.output, run_protein_id_filtered, run_peptide_ids_filtered)


if __name__ == "__main__":
    main()
