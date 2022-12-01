import os
import sys
import logging
import csv
import re
import vcf
import argparse
import urllib2
import itertools
import pandas as pd
import numpy as np
import Fred2.Core.Generator as generator
import math

from collections import defaultdict
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.Core.Variant import Variant, VariationType, MutationSyntax
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.ADBAdapter import EIdentifierTypes
from Fred2.IO.UniProtAdapter import UniProtDB
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.IO import FileReader
from Bio import SeqUtils

from datetime import datetime
from string import Template

__author__ = "mohr, walzer"
VERSION = "1.1"

ID_SYSTEM_USED = EIdentifierTypes.ENSEMBL
transcriptProteinMap = {}
transcriptSwissProtMap = {}


REPORT_TEMPLATE = """
###################################################################
        EPAA - EPITOPE PREDICTION REPORT
###################################################################
Persons in Charge: Christopher Mohr, Mathias Walzer
Date: $date
Pipeline Version: 1.1
Workflow Version: 1.1
Sample ID: $sample
Alleles
-------------
$alleles
Used Prediction Methods
-------------
$methods
Used Reference
-------------
$reference
Binding Assessment Criteria
-------------
Syfpeithi predictions: prediction score > half max score of corresponding allele
netMHC/netMHCpan predictions: affinity (as IC50 value in nM) <= 500
Additional Steps
-------------
NO filtering for peptide input.
Filtering of self-peptides (Reviewed (Swiss-Prot) UP000005640 uniprot-all.fasta.gz - 29/02/16, ENSEMBL release 84 Homo_sapiens.GRCh38.pep.all.fa.gz - 27/04/2016)
When personalized protein sequences are provided, peptides will be filtered against those as well.
Stats
-------------
Number of Variants: $variants
Number of Peptides: $peptides
Number of Peptides after Filtering: $filter
Number of Predictions: $predictions
Number of Predicted Binders: $binders
Number of Predicted Non-Binders: $nonbinders
Number of Binding Peptides: $uniquebinders
Number of Non-Binding Peptides: $uniquenonbinders
Contacts
-------------
mohr@informatik.uni-tuebingen.de
walzer@informatik.uni-tuebingen.de
University of Tuebingen, Applied Bioinformatics,
Center for Bioinformatics, Quantitative Biology Center,
and Dept. of Computer Science,
Sand 14, 72076 Tuebingen, Germany
"""


def get_fred2_annotation(vt, p, r, alt):
    if vt == VariationType.SNP:
        return p, r, alt
    elif vt == VariationType.DEL or vt == VariationType.FSDEL:
        # more than one observed ?
        if alt != "-":
            alternative = "-"
            reference = r[len(alt) :]
            position = p + len(alt)
        else:
            return p, r, alt
    elif vt == VariationType.INS or vt == VariationType.FSINS:
        if r != "-":
            position = p
            reference = "-"
            if alt != "-":
                alt_new = alt[len(r) :]
                alternative = alt_new
            else:
                alternative = str(alt)
        else:
            return p, r, alt

    return position, reference, alternative


def read_vcf(filename, pass_only=True):
    """
    reads vcf files
    returns a list of FRED2 variants
    :param filename: /path/to/file
    :return: list of FRED2 variants
    """
    global ID_SYSTEM_USED

    vl = list()
    with open(filename, "rb") as tsvfile:
        vcf_reader = vcf.Reader(tsvfile)
        vl = [r for r in vcf_reader]

    dict_vars = {}
    list_vars = []
    transcript_ids = []
    genotye_dict = {"het": False, "hom": True, "ref": True}

    for num, record in enumerate(vl):
        c = record.CHROM.strip("chr")
        p = record.POS - 1
        variation_dbid = record.ID
        r = str(record.REF)
        v_list = record.ALT
        f = record.FILTER
        if pass_only and f:
            continue

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
        """
        vt = VariationType.UNKNOWN
        if record.is_snp:
            vt = VariationType.SNP
        elif record.is_indel:
            if len(v_list) % 3 == 0:  # no frameshift
                if record.is_deletion:
                    vt = VariationType.DEL
                else:
                    vt = VariationType.INS
            else:  # frameshift
                if record.is_deletion:
                    vt = VariationType.FSDEL
                else:
                    vt = VariationType.FSINS
        gene = ""

        for alt in v_list:
            isHomozygous = False
            if "HOM" in record.INFO:
                isHomozygous = record.INFO["HOM"] == 1
            elif "SGT" in record.INFO:
                zygosity = record.INFO["SGT"].split("->")[1]
                if zygosity in genotye_dict:
                    isHomozygous = genotye_dict[zygosity]
                else:
                    if zygosity[0] == zygosity[1]:
                        isHomozygous = True
                    else:
                        isHomozygous = False
            else:
                for sample in record.samples:
                    if "GT" in sample.data:
                        isHomozygous = sample.data["GT"] == "1/1"

            if record.INFO["ANN"]:
                isSynonymous = False
                coding = dict()
                types = []
                for annraw in record.INFO["ANN"]:  # for each ANN only add a new coding! see GSvar
                    annots = annraw.split("|")
                    (
                        obs,
                        a_mut_type,
                        impact,
                        a_gene,
                        a_gene_id,
                        feature_type,
                        transcript_id,
                        exon,
                        tot_exon,
                        trans_coding,
                        prot_coding,
                        cdna,
                        cds,
                        aa,
                        distance,
                        warnings,
                    ) = annots
                    types.append(a_mut_type)

                    tpos = 0
                    ppos = 0
                    positions = ""

                    # get cds/protein positions and convert mutation syntax to FRED2 format
                    if trans_coding != "":
                        positions = re.findall(r"\d+", trans_coding)
                        ppos = int(positions[0]) - 1

                    if prot_coding != "":
                        positions = re.findall(r"\d+", prot_coding)
                        tpos = int(positions[0]) - 1

                    isSynonymous = a_mut_type == "synonymous_variant"

                    gene = a_gene_id
                    # there are no isoforms in biomart
                    transcript_id = transcript_id.split(".")[0]

                    if "NM" in transcript_id:
                        ID_SYSTEM_USED = EIdentifierTypes.REFSEQ

                    # take online coding variants into account, FRED2 cannot deal with stopgain variants right now
                    if not prot_coding or "stop_gained" in a_mut_type:
                        continue

                    coding[transcript_id] = MutationSyntax(transcript_id, ppos, tpos, trans_coding, prot_coding)
                    transcript_ids.append(transcript_id)

                if coding:
                    pos, reference, alternative = get_fred2_annotation(vt, p, r, str(alt))
                    var = Variant(
                        "line" + str(num), vt, c, pos, reference, alternative, coding, isHomozygous, isSynonymous
                    )
                    var.gene = gene
                    var.log_metadata("vardbid", variation_dbid)
                    dict_vars[var] = var
                    list_vars.append(var)

    transToVar = {}

    # fix because of memory/timing issues due to combinatoric explosion
    for v in list_vars:
        for trans_id in v.coding.iterkeys():
            transToVar.setdefault(trans_id, []).append(v)

    for tId, vs in transToVar.iteritems():
        if len(vs) > 10:
            for v in vs:
                vs_new = Variant(v.id, v.type, v.chrom, v.genomePos, v.ref, v.obs, v.coding, True, v.isSynonymous)
                vs_new.gene = v.gene
                for m in ["vardbid"]:
                    vs_new.log_metadata(m, v.get_metadata(m))
                dict_vars[v] = vs_new

    return dict_vars.values(), transcript_ids
