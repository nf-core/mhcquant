#!/usr/bin/env python

import argparse

def main():
    model = argparse.ArgumentParser(description='Neoepitope prediction for TargetInspector.')

    model.add_argument(
        '-i', '--input',
        type=str,
        help='predicted class 2 neoepitopes'
    )

    model.add_argument(
        '-n', '--neoepitopes',
        type=str,
        help='neoepitopes file'
    )

    model.add_argument(
        '-o', '--output',
        type=str,
        help='postprocessed output files'
    )

    args = model.parse_args()



if __name__ == '__main__':
    main()