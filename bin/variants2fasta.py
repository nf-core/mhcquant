#!/usr/bin/env python
import time
import sys
import argparse

from Fred2.Core import Protein, Peptide, Allele, MutationSyntax, Variant
from Fred2.Core.Variant import VariationType
from Fred2.IO import read_lines, MartsAdapter, read_annovar_exonic
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.Core import (
    generate_peptides_from_proteins,
    generate_peptides_from_variants,
    generate_transcripts_from_variants,
    generate_proteins_from_transcripts,
)
from Fred2.IO.ADBAdapter import EIdentifierTypes, EAdapterFields
from vcf_reader import read_vcf


MARTDBURL = {
    "GRCH37": "http://grch37.ensembl.org/biomart/martservice?query=",
    "GRCH38": "http://www.ensembl.org/biomart/martservice?query=",
}  # is corrently set to GRCh38


def read_variant_effect_predictor(file, gene_filter=None):
    """
    Reads a VCF (v4.1) file generatede by variant effect predictor and generates variant objects
    :param str file: Path to vcf file
    :param list gene_filter: List of proteins (in HGNC) of inerrest. Variants are filter according to this list
    :return: list(Variant) - a list of Fred2.Core.Variant objects
    """
    vars = []

    def get_type(ref, alt):
        """
        returns the variant type
        """
        if len(ref) == 1 and len(alt) == 1:
            return VariationType.SNP
        if len(ref) > 0 and len(alt) == 0:
            if len(ref) % 3 == 0:
                return VariationType.DEL
            else:
                return VariationType.FSDEL
        if len(ref) == 0 and len(alt) > 0:
            if len(alt) % 3 == 0:
                return VariationType.INS
            else:
                return VariationType.FSINS
        return VariationType.UNKNOWN

    coding_types = set(
        [
            "3_prime_UTR_variant",
            "5_prime_UTR_variant",
            "start_lost",
            "stop_gained",
            "frameshift_variant",
            "start_lost",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
            "splice_region_variant",
            "incomplete_terminal_codon_variant",
            "stop_retained_variant",
            "synonymous_variant",
            "coding_sequence_variant",
        ]
    )

    with open(file, "r") as f:
        for i, l in enumerate(f):

            # skip comments
            if l.startswith("#") or l.strip() == "":
                continue

            chrom, gene_pos, var_id, ref, alt, _, filter_flag, info = l.strip().split("\t")[:8]
            coding = {}
            isSynonymous = False

            for co in info.split(","):
                # Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|TSL|APPRIS|SIFT|PolyPhen|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">
                (
                    _,
                    var_type,
                    _,
                    gene,
                    _,
                    transcript_type,
                    transcript_id,
                    _,
                    _,
                    _,
                    _,
                    _,
                    transcript_pos,
                    _,
                    prot_pos,
                    aa_mutation,
                ) = co.strip().split("|")[:16]
                HGNC_ID = co.strip().split("|")[22]

                # pass every other feature type except Transcript (RegulatoryFeature, MotifFeature.)
                # pass genes that are uninterresting for us
                if transcript_type != "Transcript" or (HGNC_ID not in gene_filter and gene_filter):
                    continue

                # pass all intronic and other mutations that do not directly influence the protein sequence
                if any(t in coding_types for t in var_type.split("&")):
                    # generate mutation syntax

                    # positioning in Fred2 is 0-based!!!
                    if transcript_pos != "":
                        coding[transcript_id] = MutationSyntax(
                            transcript_id,
                            int(transcript_pos.split("/")[0]) - 1,
                            -1 if prot_pos == "" else int(prot_pos) - 1,
                            co,
                            "",
                            geneID=HGNC_ID,
                        )

                # is variant synonymous?
                isSynonymous = any(t == "synonymous_variant" for t in var_type.split("&"))
            if coding:
                vars.append(
                    Variant(
                        var_id,
                        get_type(ref, alt),
                        chrom,
                        int(gene_pos),
                        ref.upper(),
                        alt.upper(),
                        coding,
                        False,
                        isSynonymous,
                    )
                )
    return vars


def main():

    model = argparse.ArgumentParser(description="Neoepitope protein fasta generation from variant vcf")

    model.add_argument("-v", "--vcf", type=str, default=None, help="Path to the vcf input file")

    model.add_argument(
        "-t",
        "--type",
        type=str,
        choices=["VEP", "ANNOVAR", "SNPEFF"],
        default="VEP",
        help="Type of annotation tool used (Variant Effect Predictor, ANNOVAR exonic gene annotation, SnpEff)",
    )

    model.add_argument("-f", "--fasta_ref", type=str, default=None, help="Path to the fasta input file")

    model.add_argument(
        "-p", "--proteins", type=str, default=None, help="Path to the protein ID input file (in HGNC-ID)"
    )

    model.add_argument(
        "-r",
        "--reference",
        type=str,
        default="GRCh38",
        help="The reference genome used for varinat annotation and calling.",
    )

    model.add_argument(
        "-fINDEL", "--filterINDEL", action="store_true", help="Filter insertions and deletions (including frameshifts)"
    )

    model.add_argument("-fFS", "--filterFSINDEL", action="store_true", help="Filter frameshift INDELs")

    model.add_argument("-fSNP", "--filterSNP", action="store_true", help="Filter SNPs")

    model.add_argument("-o", "--output", type=str, required=True, help="Path to the output file")

    args = model.parse_args()

    martDB = MartsAdapter(biomart=MARTDBURL[args.reference.upper()])

    if args.vcf is None:
        sys.stderr.write("At least a vcf file or a protein id file has to be provided.\n")
        return -1

    # if vcf file is given: generate variants and filter them if HGNC IDs ar given
    if args.vcf is not None:
        protein_ids = []
        if args.proteins is not None:
            with open(args.proteins, "r") as f:
                for l in f:
                    l = l.strip()
                    if l != "":
                        protein_ids.append(l)

        if args.type == "VEP":
            variants = read_variant_effect_predictor(args.vcf, gene_filter=protein_ids)

        elif args.type == "SNPEFF":
            variants = read_vcf(args.vcf)[0]

        else:
            variants = read_annovar_exonic(args.vcf, gene_filter=protein_ids)

        if args.filterSNP:
            variants = filter(lambda x: x.type != VariationType.SNP, variants)

        if args.filterINDEL:
            variants = filter(
                lambda x: x.type
                not in [VariationType.INS, VariationType.DEL, VariationType.FSDEL, VariationType.FSINS],
                variants,
            )

        if args.filterFSINDEL:
            variants = filter(lambda x: x.type not in [VariationType.FSDEL, VariationType.FSINS], variants)

        if not variants:
            sys.stderr.write("No variants left after filtering. Please refine your filtering criteria.\n")
            return -1

        # generate transcripts
        transcripts = generate_transcripts_from_variants(variants, martDB, EIdentifierTypes.ENSEMBL)

        # generate proteins
        proteins = generate_proteins_from_transcripts(transcripts)

        # write fasta file
        with open(args.output, "w") as f:
            for p in proteins:
                f.write(">" + str(p.transcript_id) + "|" + str(p.vars) + "_var_" + "\n")
                f.write(str(p) + "\n")

        # concatenate fasta file with fasta reference
        with open(args.output) as op:
            opr1 = op.read()

        with open(args.fasta_ref) as op:
            opr2 = op.read()

        concat = opr1 + opr2

        with open(args.output, "w") as op:
            op.write(concat)

    else:
        sys.stderr.write("At least a vcf file or a protein id file has to be provided.\n")
        return -1

    return 0


if __name__ == "__main__":
    sys.exit(main())
