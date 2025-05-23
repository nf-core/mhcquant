#!/usr/bin/env python

"""
Script to convert SDRF files to MHCquant-compatible samplesheet format.
"""

import argparse
import pandas as pd
import sys
import os
from sdrf.sdrf import SdrfDataFrame

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Convert SDRF files to MHCquant-compatible samplesheet format.')
    parser.add_argument('-s', '--sdrf', required=True, help='Path to SDRF file')
    parser.add_argument('-o', '--output', required=True, help='Path to output samplesheet file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    return parser.parse_args()

def convert_sdrf_to_mhcquant(sdrf_file, output_file):
    """
    Convert SDRF file to MHCquant-compatible samplesheet format.
    
    Parameters:
    -----------
    sdrf_file : str
        Path to SDRF file
    output_file : str
        Path to output samplesheet file
    """
    try:
        # Parse SDRF file
        sdrf_df = SdrfDataFrame.parse_sdrf(sdrf_file)
        
        # Create output dataframe with required columns
        output_df = pd.DataFrame(columns=['ID', 'Sample', 'Condition', 'ReplicateFileName'])
        
        # Extract necessary information from SDRF
        for idx, row in sdrf_df.iterrows():
            # Get source name as sample
            sample = row.get('source name', f'sample_{idx}')
            
            # Get condition from factor value or characteristics
            condition = None
            for col in row.index:
                if col.startswith('factor value['):
                    condition = row[col]
                    break
            
            if condition is None:
                # If no factor value is found, use organism part or other characteristic
                for col in row.index:
                    if col.startswith('characteristics[organism part]'):
                        condition = row[col]
                        break
                
                if condition is None:
                    # Default condition if nothing else is found
                    condition = 'unknown'
            
            # Get file path
            file_path = None
            for col in row.index:
                if col == 'comment[data file]':
                    file_path = row[col]
                    break
            
            if file_path is None:
                continue
            
            # Add to output dataframe
            output_df = output_df._append({
                'ID': idx + 1,
                'Sample': sample,
                'Condition': condition,
                'ReplicateFileName': file_path
            }, ignore_index=True)
        
        # Write to output file
        output_df.to_csv(output_file, sep='\t', index=False)
        print(f"Successfully converted {sdrf_file} to {output_file}")
        
    except Exception as e:
        print(f"Error converting SDRF file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """Main function."""
    args = parse_args()
    convert_sdrf_to_mhcquant(args.sdrf, args.output)

if __name__ == '__main__':
    main()