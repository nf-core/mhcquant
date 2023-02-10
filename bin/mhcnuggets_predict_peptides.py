#!/usr/bin/env python
from mhcnuggets.src.predict import predict
import argparse
import logging
import sys

# logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("MHCNuggets Predict Peptides")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)

supported_alleles_class_2 = [
    "HLA-DPA1*01:03-DPB1*02:01",
    "HLA-DPA1*01:03-DPB1*03:01",
    "HLA-DPA1*01:03-DPB1*04:01",
    "HLA-DPA1*01:03-DPB1*04:02",
    "HLA-DPA1*02:01-DPB1*01:01",
    "HLA-DPA1*02:01-DPB1*05:01",
    "HLA-DPA1*02:02-DPB1*05:01",
    "HLA-DPA1*03:01-DPB1*04:02",
    "HLA-DPB1*01:01",
    "HLA-DPB1*02:01",
    "HLA-DPB1*03:01",
    "HLA-DPB1*04:01",
    "HLA-DPB1*04:02",
    "HLA-DPB1*05:01",
    "HLA-DPB1*09:01",
    "HLA-DPB1*11:01",
    "HLA-DPB1*14:01",
    "HLA-DPB1*20:01",
    "HLA-DQA1*01:01",
    "HLA-DQA1*01:01-DQB1*05:01",
    "HLA-DQA1*01:01-DQB1*05:03",
    "HLA-DQA1*01:02",
    "HLA-DQA1*01:02-DQB1*05:01",
    "HLA-DQA1*01:02-DQB1*05:02",
    "HLA-DQA1*01:02-DQB1*06:02",
    "HLA-DQA1*01:02-DQB1*06:04",
    "HLA-DQA1*01:03-DQB1*03:02",
    "HLA-DQA1*01:03-DQB1*06:01",
    "HLA-DQA1*01:03-DQB1*06:03",
    "HLA-DQA1*01:04-DQB1*05:03",
    "HLA-DQA1*02:01-DQB1*02:01",
    "HLA-DQA1*02:01-DQB1*02:02",
    "HLA-DQA1*02:01-DQB1*03:01",
    "HLA-DQA1*02:01-DQB1*03:03",
    "HLA-DQA1*02:01-DQB1*04:02",
    "HLA-DQA1*03:01",
    "HLA-DQA1*03:01-DQB1*02:01",
    "HLA-DQA1*03:01-DQB1*03:01",
    "HLA-DQA1*03:01-DQB1*03:02",
    "HLA-DQA1*03:01-DQB1*04:01",
    "HLA-DQA1*03:02-DQB1*03:01",
    "HLA-DQA1*03:02-DQB1*03:03",
    "HLA-DQA1*03:02-DQB1*04:01",
    "HLA-DQA1*03:03-DQB1*04:02",
    "HLA-DQA1*04:01-DQB1*04:02",
    "HLA-DQA1*05:01",
    "HLA-DQA1*05:01-DQB1*02:01",
    "HLA-DQA1*05:01-DQB1*03:01",
    "HLA-DQA1*05:01-DQB1*03:02",
    "HLA-DQA1*05:01-DQB1*03:03",
    "HLA-DQA1*05:01-DQB1*04:02",
    "HLA-DQA1*05:05-DQB1*03:01",
    "HLA-DQA1*06:01-DQB1*04:02",
    "HLA-DQB1*02:01",
    "HLA-DQB1*02:02",
    "HLA-DQB1*03:01",
    "HLA-DQB1*03:02",
    "HLA-DQB1*03:19",
    "HLA-DQB1*04:02",
    "HLA-DQB1*05:01",
    "HLA-DQB1*05:02",
    "HLA-DQB1*05:03",
    "HLA-DQB1*06:02",
    "HLA-DQB1*06:03",
    "HLA-DQB1*06:04",
    "HLA-DRA0*10:1-DRB1*01:01",
    "HLA-DRA0*10:1-DRB1*03:01",
    "HLA-DRA0*10:1-DRB1*04:01",
    "HLA-DRA0*10:1-DRB1*04:04",
    "HLA-DRA0*10:1-DRB1*07:01",
    "HLA-DRA0*10:1-DRB1*08:01",
    "HLA-DRA0*10:1-DRB1*09:01",
    "HLA-DRA0*10:1-DRB1*11:01",
    "HLA-DRA0*10:1-DRB1*13:01",
    "HLA-DRA0*10:1-DRB1*14:54",
    "HLA-DRA0*10:1-DRB1*15:01",
    "HLA-DRA0*10:1-DRB3*01:01",
    "HLA-DRA0*10:1-DRB3*02:02",
    "HLA-DRA0*10:1-DRB3*03:01",
    "HLA-DRA0*10:1-DRB4*01:03",
    "HLA-DRA0*10:1-DRB5*01:01",
    "HLA-DRB1*01:01",
    "HLA-DRB1*01:02",
    "HLA-DRB1*01:03",
    "HLA-DRB1*03:01",
    "HLA-DRB1*03:02",
    "HLA-DRB1*03:03",
    "HLA-DRB1*03:04",
    "HLA-DRB1*03:05",
    "HLA-DRB1*04:01",
    "HLA-DRB1*04:02",
    "HLA-DRB1*04:03",
    "HLA-DRB1*04:04",
    "HLA-DRB1*04:05",
    "HLA-DRB1*04:06",
    "HLA-DRB1*04:07",
    "HLA-DRB1*04:11",
    "HLA-DRB1*07:01",
    "HLA-DRB1*08:01",
    "HLA-DRB1*08:02",
    "HLA-DRB1*08:03",
    "HLA-DRB1*08:04",
    "HLA-DRB1*09:01",
    "HLA-DRB1*10:01",
    "HLA-DRB1*11:01",
    "HLA-DRB1*11:02",
    "HLA-DRB1*11:03",
    "HLA-DRB1*11:04",
    "HLA-DRB1*12:01",
    "HLA-DRB1*12:02",
    "HLA-DRB1*13:01",
    "HLA-DRB1*13:02",
    "HLA-DRB1*13:03",
    "HLA-DRB1*13:04",
    "HLA-DRB1*13:05",
    "HLA-DRB1*14:01",
    "HLA-DRB1*14:02",
    "HLA-DRB1*15:01",
    "HLA-DRB1*15:02",
    "HLA-DRB1*15:03",
    "HLA-DRB1*16:01",
    "HLA-DRB1*16:02",
    "HLA-DRB3*01:01",
    "HLA-DRB3*02:02",
    "HLA-DRB3*03:01",
    "HLA-DRB4*01:01",
    "HLA-DRB4*01:03",
    "HLA-DRB5*01:01",
    "HLA-DRB5*01:02",
]

flatten = lambda l: [item for sublist in l for item in sublist]


def convert_alleles_mhcnuggets_format(alleles):
    return [allele.replace("*", "") for allele in alleles]


def parse_alleles(allele_input):
    alleles = allele_input.split(";")
    supp_alleles = convert_alleles_mhcnuggets_format(list(set(alleles).intersection(supported_alleles_class_2)))

    return supp_alleles


def main():
    model = argparse.ArgumentParser(description="MHCNuggets binding prediction")

    model.add_argument("-p", "--peptides", type=str, help="mhcnuggets input")

    model.add_argument("-a", "--alleles", type=str, help="class 2 alleles")

    model.add_argument("-o", "--output", type=str, help="mhcnuggets output")

    args = model.parse_args()

    if open(args.peptides).readlines() != []:
        supp_alleles = parse_alleles(args.alleles)

        for allele in supp_alleles:
            predict(class_="II", peptides_path=args.peptides, mhc=allele, output=allele + args.output)

    else:
        op = open("predicted_neoepitopes_class_2", "w")
        op.close()


if __name__ == "__main__":
    main()
