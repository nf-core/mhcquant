#!/usr/bin/env python

# Written by Jana Hoffmann under the MIT license

import yaml
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--seq_column', type=str, help='Column containing peptides.')
    parser.add_argument('--protacc_column', type=str, help='Column containing protein accession.')
    parser.add_argument('--intensity_column', type=str, help='Column containing peptide intensity.')
    parser.add_argument('--start_column', type=str, help='Column containing peptide start position.')
    parser.add_argument('--end_column', type=str, help='Column containing peptide end position.')
    parser.add_argument('--mod_pattern', type=str, help='Delimiter separating peptide sequence from modifications.')
    parser.add_argument('--delimiter', type=str, help='Delimiter separating multiple entries in one cell.')
    parser.add_argument('--min_overlap', type=int, help='Minimal overlap of two peptides to be grouped together.')
    parser.add_argument('--max_step_size', type=int, help='Maximal difference in the start position of two peptides to be still grouped together.')
    parser.add_argument('--min_epi_length', type=int, help='Minimal epitope length.')
    parser.add_argument('--prefix', type=str, help='The prefix for the params.yml.')

    return parser.parse_args()

def create_yml(prefix, seq_column, protacc_column, intensity_column, start_column, end_column, mod_pattern, delimiter, min_overlap, max_step_size, min_epi_length):
    """Create a yml file from all input parameters of epicore.

    Args:
        prefix (str):           The prefix for the params.yml
        seq_column (str):       The column containing the peptide sequences
        protacc_column (str):   The column containing the protein accessions
        intensity_column (str): Column containing peptide intensity
        start_column (str):     Column containing peptide start position
        end_column (str):       Column containing peptide end position
        mod_pattern (str):      Delimiter separating peptide sequence from modifications
        delimiter (str):        Delimiter separating multiple entries in one cell
        min_overlap (str):      Minimal overlap of two peptides to be grouped together
        max_step_size (str):    Maximal difference in the start position of two peptides to be still grouped together
        min_epi_length (str):   Minimal epitope length

    """
    yaml_list = {'parameters':{'seq_column':seq_column, 'protacc_column':protacc_column, 'intensity_column':intensity_column, 'start_column': start_column, 'end_column':end_column, 'mod_pattern':mod_pattern, 'delimiter': delimiter, 'min_overlap': min_overlap, 'max_step_size': max_step_size, 'min_epi_length': min_epi_length, 'report': '', 'out_dir': '.', 'prot_accession': ''}}

    with open(f'{prefix}_epicore_params.yml', 'w') as yaml_file:
        yaml.dump(yaml_list, yaml_file, default_flow_style=False)

def main():
    args = parse_arguments()
    seq_column = args.seq_column
    protacc_column = args.protacc_column
    intensity_column = args.intensity_column
    start_column = args.start_column
    end_column = args.end_column
    mod_pattern = args.mod_pattern
    delimiter = args.delimiter
    min_overlap = args.min_overlap
    max_step_size = args.max_step_size
    min_epi_length = args.min_epi_length
    prefix = args.prefix

    create_yml(prefix, seq_column, protacc_column, intensity_column, start_column, end_column, mod_pattern, delimiter, min_overlap, max_step_size, min_epi_length)

if __name__ == '__main__':
    main()