#!/usr/bin/env python
from mhcflurry import Class1AffinityPredictor
import logging
import sys

# logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
LOG = logging.getLogger("MHCFlurry Predict mztab for filtering")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)

# List of alleles supported by mhclurry
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
print(sys.argv[-4])
alleles = sys.argv[-4].split(";")
print(alleles)

# extract and verify alleles
unsupported_alleles = [a for a in alleles if a not in supported_alleles]
alleles = [a for a in alleles if a in supported_alleles]

if unsupported_alleles:
    for allele in unsupported_alleles:
        LOG.warning("Allele: " + allele + " is not supported by MHCFlurry!")
if not alleles:
    LOG.warning("Submitted alleles are not supported or formatting of input.tsv is not correct!")


# read identified peptides with q-value < threshold
mztab = open(sys.argv[-3])
mztab_read = mztab.readlines()
mztab.close()
seqs = [l.split()[1] for l in mztab_read if l.startswith("PSM")]
seqs_new_smaller_qval = list(set(seqs))

# read all remaining peptides with q-value > threshold
mztab = open(sys.argv[-2])
mztab_read = mztab.readlines()
mztab.close()
seqs = [l.split()[1] for l in mztab_read if l.startswith("PSM") if l.split()[1] not in seqs_new_smaller_qval]
seqs_new_greater_qval = list(set(seqs))
seqs_new_greater_qval = [
    s
    for s in seqs_new_greater_qval
    if 7 < len(s) < 13 and not "U" in s and not "X" in s and not "Z" in s and not "J" in s and not "B" in s
]

# call mhcflurry
# subprocess.call("mhcflurry-predict --peptides {p} --alleles {a} --out {o}".format(p=" ".join(seqs_new), a=" ".join(alleles), o=sys.argv[-1]))
seqs_filtered = []
for allele in alleles:
    print(allele)
    predictor = Class1AffinityPredictor.load()
    df_pred = predictor.predict_to_dataframe(allele=allele, peptides=seqs_new_greater_qval)
    seqs_filtered += df_pred[df_pred["prediction"] <= float(sys.argv[-5])]["peptide"].values.tolist()

# merge sequence lists and append decoys
seqs_new_all = list(set(seqs_new_smaller_qval + seqs_filtered))
seqs_new_all = seqs_new_all + [s[::-1] for s in seqs_new_all]

# write idXML for filtering
op = open(sys.argv[-1], "w")
op.write('<PeptideIdentification score_type="q-value" higher_score_better="false">' + "\n")
for pep in seqs_new_all:
    op.write("\t\t\t" + '<PeptideHit sequence="' + pep + '" score="0" charge="0" >' + "\n")
    op.write("\t\t\t" + "</PeptideHit>" + "\n")
op.write("</PeptideIdentification>" + "\n")
op.close()
