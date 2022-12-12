#!/usr/bin/env python
import pandas as pd
import numpy as np
import logging
import sys

from mhcflurry import Class1AffinityPredictor

# logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("Neoepitope Binding Prediction")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)

# List of alleles supported by mhcflurry
supported_alleles = (
    "A*01:01,A*02:01,A*02:02,A*02:03,A*02:05,A*02:06,A*02:07,A*02:11,A*02:12,A*02:16,A*02:17,A*02:19,"
    "A*02:50,A*03:01,A*11:01,A*23:01,A*24:02,A*24:03,A*25:01,A*26:01,A*26:02,A*26:03,A*29:02,A*30:01,"
    "A*30:02,A*31:01,A*32:01,A*33:01,A*66:01,A*68:01,A*68:02,A*68:23,A*69:01,A*80:01,B*07:01,B*07:02,"
    "B*08:01,B*08:02,B*08:03,B*14:02,B*15:01,B*15:02,B*15:03,B*15:09,B*15:17,B*18:01,B*27:02,B*27:03,"
    "B*27:04,B*27:05,B*27:06,B*35:01,B*35:03,B*37:01,B*38:01,B*39:01,B*39:06,B*40:01,B*40:02,B*42:01,"
    "B*44:02,B*44:03,B*45:01,B*46:01,B*48:01,B*51:01,B*53:01,B*54:01,B*57:01,B*58:01,B*83:01,C*03:03,"
    "C*04:01,C*05:01,C*06:02,C*07:02,C*08:02,C*12:03,C*14:02,C*15:02".split(",")
)

# read provided allotypes
alleles = sys.argv[-3].split(";")

# extract and verify alleles
unsupported_alleles = [a for a in alleles if a not in supported_alleles]
alleles = [a for a in alleles if a in supported_alleles]

if unsupported_alleles:
    for allele in unsupported_alleles:
        LOG.warning("Allele: " + allele + " is not supported by MHCFlurry!")
if not alleles:
    LOG.warning("Submitted alleles are not supported or formatting of input.tsv is not correct!")

flatten = lambda l: [item for sublist in l for item in sublist]
# read identified peptides
neoepitopes = [line.rstrip("\n").strip().split(",") for line in open(sys.argv[-2])][1:]
neoepitopes = flatten(neoepitopes)
seqs_to_geneID = dict(zip(neoepitopes[::2], neoepitopes[1::2]))

if len(seqs_to_geneID) > 0:
    # call mhcflurry
    for allele in alleles:
        predictor = Class1AffinityPredictor.load()
        df_pred = predictor.predict_to_dataframe(allele=allele, peptides=seqs_to_geneID.keys())
        df_pred.insert(1, "geneID", pd.Series(np.array(seqs_to_geneID.values())))
        df_pred.to_csv(allele + "_" + sys.argv[-1])
else:
    op = open(sys.argv[-1], "w")
    op.close()
