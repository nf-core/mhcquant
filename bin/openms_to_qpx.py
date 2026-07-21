#!/usr/bin/env python
# Written by Jonas Scheid and released under MIT license.

"""Convert mhcquant consensusXML file(s) into qpx psm.parquet and feature.parquet."""

import argparse
import logging
import os
import re

import pandas as pd
import pyopenms as oms
from qpx.writers import FeatureWriter, PsmWriter

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

RAW_FILE_EXTENSIONS = (".mzML", ".mzml", ".raw", ".RAW")
PSM_DEDUP_KEYS = ["sequence", "charge", "run_file_name"]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert mhcquant consensusXML file(s) into qpx psm.parquet and feature.parquet."
    )
    parser.add_argument(
        "--consensusxml",
        required=True,
        nargs="+",
        help="One or more consensusXML files (e.g. one per Sample+Condition group) to convert.",
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Project accession, used as the output file name prefix.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for <accession>.psm.parquet and <accession>.feature.parquet.",
    )
    return parser.parse_args()


def normalize_run_name(filename: str) -> str:
    """qpx run_file_name is the raw data file stem; strip mhcquant's _cleaned/_aligned preprocessing suffixes to match the SDRF."""
    name = filename.rsplit("/", 1)[-1]
    for ext in RAW_FILE_EXTENSIONS:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    name = re.sub(r"_aligned$", "", name)
    name = re.sub(r"_cleaned$", "", name)
    return name


def to_proforma(aa_sequence: "oms.AASequence") -> str:
    """Render an OpenMS AASequence (incl. N-/C-terminal mods) as ProForma using UniMod accessions."""
    proforma = ""
    n_term = aa_sequence.getNTerminalModification()
    if n_term and n_term.getUniModAccession():
        proforma += "[" + n_term.getUniModAccession().replace("UniMod:", "UNIMOD:") + "]-"
    for i in range(aa_sequence.size()):
        residue = aa_sequence.getResidue(i)
        proforma += residue.getOneLetterCode()
        mod = residue.getModification()
        if mod and mod.getUniModAccession():
            proforma += "[" + mod.getUniModAccession().replace("UniMod:", "UNIMOD:") + "]"
    c_term = aa_sequence.getCTerminalModification()
    if c_term and c_term.getUniModAccession():
        proforma += "-[" + c_term.getUniModAccession().replace("UniMod:", "UNIMOD:") + "]"
    return proforma


def scan_numbers(peptide_id: "oms.PeptideIdentification") -> list:
    spectrum_ref = (
        peptide_id.getMetaValue("spectrum_reference") if peptide_id.metaValueExists("spectrum_reference") else ""
    )
    match = re.search(r"scan=(\d+)", spectrum_ref or "")
    return [int(match.group(1))] if match else [0]


def rows_from_consensusxml(path: str) -> tuple:
    """Extract (psm_rows, feature_rows): one psm row per identified consensus feature, one feature row per contributing run."""
    cm = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(path, cm)
    run_by_map_index = {idx: normalize_run_name(header.filename) for idx, header in cm.getColumnHeaders().items()}
    if any(not name for name in run_by_map_index.values()):
        # mhcquant single-run consensus leaves the column header filename empty; the run name is in primaryMSRunPath
        primary_runs = []
        for protein_id in cm.getProteinIdentifications():
            paths = []
            protein_id.getPrimaryMSRunPath(paths)
            primary_runs = [normalize_run_name(p.decode() if isinstance(p, bytes) else p) for p in paths]
            if primary_runs:
                break
        for idx in sorted(run_by_map_index):
            if not run_by_map_index[idx]:
                run_by_map_index[idx] = primary_runs[idx] if idx < len(primary_runs) else (primary_runs[0] if primary_runs else "")

    psm_rows, feature_rows = [], []
    for cf in cm:
        peptide_ids = cf.getPeptideIdentifications()
        if not peptide_ids or not peptide_ids[0].getHits():
            continue
        pid = peptide_ids[0]
        hit = pid.getHits()[0]
        sequence = hit.getSequence()
        charge = hit.getCharge()
        target_decoy = hit.getMetaValue("target_decoy") if hit.metaValueExists("target_decoy") else "target"
        protein_accessions = [ev.getProteinAccession() for ev in hit.getPeptideEvidences()]
        anchor_protein = protein_accessions[0] if protein_accessions else "unknown"

        feature_handles = cf.getFeatureList()
        id_run = run_by_map_index.get(feature_handles[0].getMapIndex(), "NA") if feature_handles else "NA"

        additional_scores = [
            {"score_name": pid.getScoreType(), "score_value": float(hit.getScore()), "higher_better": pid.isHigherScoreBetter()}
        ]
        if hit.metaValueExists("q-value"):
            additional_scores.append(
                {"score_name": "q-value", "score_value": float(hit.getMetaValue("q-value")), "higher_better": False}
            )

        base = dict(
            sequence=sequence.toUnmodifiedString(),
            peptidoform=to_proforma(sequence),
            charge=charge,
            is_decoy=str(target_decoy).startswith("decoy"),
            calculated_mz=float(sequence.getMZ(charge)),
            observed_mz=float(pid.getMZ()),
            scan=scan_numbers(pid),
            rt=float(pid.getRT()),
            additional_scores=additional_scores,
        )
        if hit.metaValueExists("MS:1001493"):
            base["posterior_error_probability"] = float(hit.getMetaValue("MS:1001493"))

        psm_rows.append({**base, "run_file_name": id_run, "protein_accessions": protein_accessions or [anchor_protein]})

        for handle in feature_handles:
            feature_rows.append(
                {
                    **base,
                    "run_file_name": run_by_map_index.get(handle.getMapIndex(), "NA"),
                    "anchor_protein": anchor_protein,
                    "id_run_file_name": id_run,
                    "intensities": [{"label": "label free sample", "intensity": float(handle.getIntensity())}],
                }
            )

    return psm_rows, feature_rows


def dedup_psm(rows: list) -> pd.DataFrame:
    """ConsensusMap iteration can revisit the same PSM via nearby co-located features; scan is a list so extract a scalar key to dedup on."""
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    scan_key = df["scan"].map(lambda scans: scans[0] if scans else -1)
    return df.assign(_scan=scan_key).drop_duplicates(subset=[*PSM_DEDUP_KEYS, "_scan"]).drop(columns="_scan")


def reindex_and_write(df: pd.DataFrame, writer) -> int:
    """qpx writers do pa.Table.from_pandas(df, schema=...): every schema column must be present, so fill missing optional ones with None."""
    for column in writer.arrow_schema.names:
        if column not in df.columns:
            df[column] = None
    writer.write_dataframe(df[writer.arrow_schema.names])
    writer.close()
    return len(df)


def main():
    args = parse_arguments()
    os.makedirs(args.outdir, exist_ok=True)

    psm_frames, feature_rows = [], []
    for path in args.consensusxml:
        logger.info("Reading %s", path)
        psm_rows, feat_rows = rows_from_consensusxml(path)
        psm_frames.append(dedup_psm(psm_rows))
        feature_rows.extend(feat_rows)

    psm_df = pd.concat(psm_frames, ignore_index=True) if psm_frames else pd.DataFrame()
    feature_df = pd.DataFrame(feature_rows)

    if psm_df.empty:
        logger.warning("No identified PSMs in any input consensusXML; writing empty psm/feature parquet.")

    n_psm = reindex_and_write(psm_df, PsmWriter(os.path.join(args.outdir, f"{args.accession}.psm.parquet")))
    n_feature = reindex_and_write(
        feature_df, FeatureWriter(os.path.join(args.outdir, f"{args.accession}.feature.parquet"))
    )
    logger.info("Wrote %d psm rows and %d feature rows to %s", n_psm, n_feature, args.outdir)


if __name__ == "__main__":
    main()
