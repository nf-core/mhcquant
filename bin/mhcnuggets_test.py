#!/usr/bin/env python
# importing the predict module
from mhcnuggets.src.predict import predict
import argparse

def main():
    model = argparse.ArgumentParser(description='Neoepitope prediction for TargetInspector.')

    model.add_argument(
        '-i', '--input',
        type=str,
        help='mhcnuggets input'
    )

    model.add_argument(
        '-o', '--output',
        type=str,
        help='mhcnuggets output'
    )
    
    args = model.parse_args()


    # similarly doing the same prediction for MHC class_II allele HLA-DRB1*01:01
    predict(class_='II',
            peptides_path=args.input, 
            mhc='HLA-DRB101:01',
            output=args.output)



if __name__ == '__main__':
    main()