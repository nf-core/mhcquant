#!/usr/bin/env python
from mhcflurry import Class1AffinityPredictor
import sys

#List of alleles supported by mhclurry
supported_alleles="A*01:01,A*02:01,A*02:02,A*02:03,A*02:05,A*02:06,A*02:07,A*02:11,A*02:12,A*02:16,A*02:17,A*02:19,A*02:50,A*03:01,A*11:01,A*23:01,A*24:02,A*24:03,A*25:01,A*26:01,A*26:02,A*26:03,A*29:02,A*30:01,A*30:02,A*31:01,A*32:01,A*33:01,A*66:01,A*68:01,A*68:02,A*68:23,A*69:01,A*80:01,B*07:01,B*07:02,B*08:01,B*08:02,B*08:03,B*14:02,B*15:01,B*15:02,B*15:03,B*15:09,B*15:17,B*18:01,B*27:02,B*27:03,B*27:04,B*27:05,B*27:06,B*35:01,B*35:03,B*37:01,B*38:01,B*39:01,B*39:06,B*40:01,B*40:02,B*42:01,B*44:02,B*44:03,B*45:01,B*46:01,B*48:01,B*51:01,B*53:01,B*54:01,B*57:01,B*58:01,B*83:01,C*03:03,C*04:01,C*05:01,C*06:02,C*07:02,C*08:02,C*12:03,C*14:02,C*15:02".split(",")

#read provided allotypes
op=open(sys.argv[-3])
alleles=op.read().split("\n")
op.close()

alleles=[a for a in alleles if a in supported_alleles]
if alleles==[]:
   print "submitted alleles are not supported or formatting of input.tsv is not correct!"

#read identified peptides
mztab=open(sys.argv[-2])
mztab_read=mztab.readlines()
mztab.close()
seqs=[l.split()[1] for l in mztab_read if l.startswith("PEP") and not 'U' in l and not 'X' in l and not 'Z' in l and not 'J' in l and not 'B' in l]
seqs_new=list(set(seqs))

#call mhcflurry
for allele in alleles:
   predictor = Class1AffinityPredictor.load()
   df_pred=predictor.predict_to_dataframe(allele=allele, peptides=seqs_new)
   df_pred.to_csv(allele + '_' + sys.argv[-1])

