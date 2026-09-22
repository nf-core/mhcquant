#!/usr/bin/env python
# Written by Jonas Scheid under the MIT license

import sys
import click
import importlib.resources
import json
import logging
from pathlib import Path
from typing import List

from ms2rescore import rescore, package_data
from ms2rescore.exceptions import RescoringError
from psm_utils.io.idxml import IdXMLReader, IdXMLWriter
from psm_utils import PSMList
import pyopenms as oms

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# MS²Rescore >=4 always rescores with ristretto. Percolator is run downstream by the pipeline on
# the feature-annotated idXML, so in that mode the ristretto scores are only written as metavalues.
RESCORING_ENGINES = ("percolator", "ristretto")
# Pipeline FDR levels mapped onto ristretto rollups. Percolator's "peptide" level is the modified
# sequence, which corresponds to ristretto's peptidoform rollup. Protein-level FDR is not exposed by
# MS²Rescore's PSM metadata, so it falls back to the peptidoform level.
FDR_LEVELS = {"psm_level_fdrs": None, "peptide_level_fdrs": "peptidoform", "protein_level_fdrs": "peptidoform"}
FEATURE_GENERATORS = ("basic", "ms2pip", "deeplc", "ms2", "im2deep")
# q-value/PEP are the idXML keys psm_utils writes for psm.qvalue/psm.pep
LEAKING_METAVALUES = {"q-value", "PEP"}
RISTRETTO_METAVALUES = {
    "ristretto_score",
    "ristretto_psm_qvalue",
    "ristretto_psm_pep",
    "ristretto_peptidoform_qvalue",
    "ristretto_peptidoform_pep",
    "ristretto_peptide_qvalue",
    "ristretto_peptide_pep",
}
LEAKING_METAVALUES |= RISTRETTO_METAVALUES
# CLI options that are consumed while building nested config sections and must not be copied
# verbatim into the top-level MS²Rescore config.
NESTED_OPTIONS = {
    "ms2pip_model",
    "ms2pip_model_dir",
    "ms2_tolerance",
    "calibration_set_size",
    "train_fdr",
    "rescoring_engine",
    "fdr_level",
    "require_precomputed_features",
}


def parse_cli_arguments_to_config(**kwargs):
    """Update default MS²Rescore config with CLI arguments"""
    config = json.load(importlib.resources.open_text(package_data, "config_default.json"))
    ms2rescore_config = config["ms2rescore"]

    for key, value in kwargs.items():
        if key in NESTED_OPTIONS:
            continue

        elif key == "feature_generators":
            feature_generators = [fgen.strip() for fgen in value.split(",") if fgen.strip()]
            unknown = set(feature_generators) - set(FEATURE_GENERATORS)
            if unknown:
                raise click.BadParameter(
                    f"Unknown feature generator(s) {sorted(unknown)}. Choose from {', '.join(FEATURE_GENERATORS)}.",
                    param_hint="--feature_generators",
                )
            # Reset feature generator dict since there might be default generators we don't want
            ms2rescore_config["feature_generators"] = {}
            if "basic" in feature_generators:
                ms2rescore_config["feature_generators"]["basic"] = {}
            if "ms2pip" in feature_generators:
                ms2rescore_config["feature_generators"]["ms2pip"] = {
                    "model": kwargs["ms2pip_model"],
                    "model_dir": kwargs["ms2pip_model_dir"],
                }
            if "deeplc" in feature_generators:
                ms2rescore_config["feature_generators"]["deeplc"] = {
                    "finetune": False,
                    "calibration_set_size": kwargs["calibration_set_size"],
                }
            if "ms2" in feature_generators:
                ms2rescore_config["feature_generators"]["ms2"] = {}
            if "im2deep" in feature_generators:
                ms2rescore_config["feature_generators"]["im2deep"] = {}

        elif key == "output_path":
            # MS²Rescore derives all its side outputs (feature names, report, tables) from this stem
            ms2rescore_config[key] = str(Path(value).with_suffix(""))

        else:
            ms2rescore_config[key] = value

    # Fragment mass tolerance is a global setting in MS²Rescore >=4 (shared by all generators)
    ms2rescore_config["tolerance_value"] = kwargs["ms2_tolerance"]
    ms2rescore_config["tolerance_mode"] = "Da"

    # Ristretto (the only rescoring engine in MS²Rescore >=4)
    ms2rescore_config["rescoring"] = {"train_fdr": kwargs["train_fdr"], "model": "svm"}
    if kwargs["require_precomputed_features"]:
        # MS²Rescore 4.0.2 mis-assigns skipped generators' features to the "psm_file" group as well, which
        # makes the HTML report fail on duplicate columns. The per-group reports already exist, so skip it.
        ms2rescore_config["write_report"] = False
    if kwargs["rescoring_engine"] == "percolator":
        logging.info(
            "Percolator rescoring engine has been specified. Ristretto q-values/PEPs are written as "
            "metavalues only; the idXML containing rescoring features is passed on to Percolator in a separate step."
        )

    return config


def rescore_idxml(
    input_file, output_file, config, rescoring_engine: str, fdr_level: str, require_precomputed_features: bool
) -> None:
    """Rescore PSMs in an idXML file and keep other information unchanged."""
    # Read PSMs
    reader = IdXMLReader(input_file)
    psm_list = reader.read_file()

    # psm_utils treats every numeric metavalue as a rescoring feature. Scores written by a previous
    # MS²Rescore pass (e.g. per-group output merged for global FDR) must not leak back in as features.
    for psm in psm_list:
        for key in LEAKING_METAVALUES & set(psm.rescoring_features):
            del psm.rescoring_features[key]

    # Ristretto overwrites the search engine score in place. Remember it so the idXML keeps the
    # original main score (and its orientation) while ristretto q-value/PEP become metavalues.
    original_scores = {id(psm): psm.score for psm in psm_list}
    features_before = {id(psm): set(psm.rescoring_features) for psm in psm_list}

    # Rescore
    try:
        rescore(config, psm_list)
        rescored = True
    except RescoringError:
        # Ristretto needs targets passing train_fdr to train. Features are already attached to the
        # PSMs at this point, so with Percolator downstream we can still hand over the features.
        rescored = False
        if rescoring_engine == "percolator":
            logging.warning(
                "Ristretto could not be trained on this input (too few confident targets). "
                "Writing MS²Rescore features without ristretto scores; Percolator will rescore downstream."
            )
        else:
            # Mirror a group in which nothing passes FDR: q-value 1 for every PSM lets the pipeline's
            # empty-group handling report and skip it instead of aborting the whole run.
            logging.warning(
                "Ristretto could not be trained on this input (too few confident targets). "
                "All PSMs are written with q-value 1, so this group will be reported as empty after FDR filtering."
            )
            for psm in psm_list:
                psm.qvalue = 1.0
                psm.pep = 1.0

    if require_precomputed_features:
        generated = [
            sorted(set(psm.rescoring_features) - features_before[id(psm)]) for psm in psm_list
        ]
        newly_generated = sorted({f for fs in generated for f in fs} - RISTRETTO_METAVALUES)
        if newly_generated:
            raise click.ClickException(
                "Expected all rescoring features to be present in the input idXML (e.g. merged per-group "
                f"MS²Rescore output for global FDR), but generators added {newly_generated[:5]}... "
                "Check that --feature_generators matches the per-group run."
            )

    for psm in psm_list:
        if rescored and psm.pep is not None:
            psm.rescoring_features["ristretto_score"] = float(psm.score)
        psm.score = original_scores[id(psm)]

    if rescored:
        apply_fdr_level(psm_list, fdr_level)
        if rescoring_engine == "ristretto" and FDR_LEVELS[fdr_level] is not None:
            collapse_to_best_psm(psm_list)

    # Keep only PSMs that were processed by all feature generators (and survived ristretto)
    peptide_ids_filtered = filter_out_artifact_psms(psm_list, reader.peptide_ids, require_pep=rescored)

    # Write
    writer = IdXMLWriter(output_file, protein_ids=reader.protein_ids, peptide_ids=peptide_ids_filtered)
    writer.write_file(psm_list)


def collapse_to_best_psm(psm_list: PSMList) -> None:
    """Keep the rollup q-value/PEP only on the best-scoring PSM per peptidoform (charge-independent).

    Mirrors PercolatorAdapter's peptide-level output, where all but the best PSM of a peptide are set
    to q-value 1 and removed by the downstream IDFilter. PSM-level ristretto values stay available as
    `ristretto_psm_*` metavalues, and quantification re-expands PSMs from the pre-filter file.
    """
    best = {}
    for psm in psm_list:
        if psm.pep is None:
            continue
        key = psm.peptidoform.modified_sequence
        score = psm.rescoring_features["ristretto_score"]
        if key not in best or score > best[key][0]:
            best[key] = (score, id(psm))
    for psm in psm_list:
        if psm.pep is None:
            continue
        if best[psm.peptidoform.modified_sequence][1] != id(psm):
            psm.qvalue = 1.0
            psm.pep = 1.0


def apply_fdr_level(psm_list: PSMList, fdr_level: str) -> None:
    """Expose ristretto's rollup q-values/PEPs and select the level written as `q-value`/`PEP`."""
    rollup = FDR_LEVELS[fdr_level]
    if fdr_level == "protein_level_fdrs":
        logging.warning("Protein-level FDR is not available from ristretto; using peptidoform-level q-values instead.")
    for psm in psm_list:
        if psm.pep is None:  # dropped by ristretto
            continue
        # Keep all levels as metavalues (written alongside the rescoring features)
        for level in ("peptidoform", "peptide"):
            for kind in ("qvalue", "pep"):
                value = psm.metadata.get(f"{level}_{kind}")
                if value is not None:
                    psm.rescoring_features[f"ristretto_{level}_{kind}"] = float(value)
        psm.rescoring_features["ristretto_psm_qvalue"] = float(psm.qvalue)
        psm.rescoring_features["ristretto_psm_pep"] = float(psm.pep)
        if rollup is not None:
            psm.qvalue = float(psm.metadata[f"{rollup}_qvalue"])
            psm.pep = float(psm.metadata[f"{rollup}_pep"])


def filter_out_artifact_psms(
    psm_list: PSMList, peptide_ids: List[oms.PeptideIdentification], require_pep: bool = True
) -> List[oms.PeptideIdentification]:
    """Filter out PeptideHits that could not be processed by all feature generators or were dropped by ristretto"""
    num_mandatory_features = max([len(psm.rescoring_features) for psm in psm_list])
    # PEP is only set by ristretto, so PSMs without it were dropped during rescoring (e.g. rank filter)
    new_psm_list = PSMList(
        psm_list=[
            psm
            for psm in psm_list
            if len(psm.rescoring_features) == num_mandatory_features and (psm.pep is not None or not require_pep)
        ]
    )

    # get differing peptidoforms of both psm lists
    psm_list_peptides = set([next(iter(psm.provenance_data.items()))[1] for psm in psm_list])
    new_psm_list_peptides = set([next(iter(psm.provenance_data.items()))[1] for psm in new_psm_list])
    not_supported_peptides = psm_list_peptides - new_psm_list_peptides

    # no need to filter if all peptides are supported
    if len(not_supported_peptides) == 0:
        return peptide_ids
    # Create new peptide ids and filter out not supported peptides
    new_peptide_ids = []
    for peptide_id in peptide_ids:
        new_hits = []
        for hit in peptide_id.getHits():
            if hit.getSequence().toString() in not_supported_peptides:
                continue
            new_hits.append(hit)
        if len(new_hits) == 0:
            continue
        peptide_id.setHits(new_hits)
        new_peptide_ids.append(peptide_id)
    logging.info(
        f"Removed {len(psm_list_peptides) - len(new_psm_list_peptides)} PSMs. Peptides not supported: {not_supported_peptides}"
    )
    return new_peptide_ids


@click.command()
@click.option(
    "-p", "--psm_file", help="Path to PSM file (PIN, mzIdentML, MaxQuant msms, X!Tandem XML, idXML)", required=True
)
@click.option(
    "-s",
    "--spectrum_path",
    help="Path to MGF/mzML spectrum file or directory with spectrum files (default: derived from identification file)",
    required=True,
)
@click.option(
    "-o", "--output_path", help="Path and stem for output file names (default: derive from identification file)"
)
@click.option("-l", "--log_level", help="Logging level (default: `info`)", default="info")
@click.option("-n", "--processes", help="Number of parallel processes available to MS²Rescore", type=int, default=16)
@click.option(
    "-fg",
    "--feature_generators",
    help=f"Comma-separated list of feature generators to use (default: `ms2pip,deeplc`). Choose from {','.join(FEATURE_GENERATORS)}",
    default="ms2pip,deeplc",
)
@click.option("-pipm", "--ms2pip_model", help="MS²PIP model (default: `Immuno-HCD`)", type=str, default="Immuno-HCD")
@click.option("-pipmdir", "--ms2pip_model_dir", help="Path to directory, which holds pre-downloaded MS²PIP models", type=str, default=None)
@click.option(
    "-ms2tol", "--ms2_tolerance", help="Fragment mass tolerance [Da](default: `0.02`)", type=float, default=0.02
)
@click.option(
    "-cs",
    "--calibration_set_size",
    help="Percentage of number of calibration set for DeepLC (default: `0.15`)",
    default=0.15,
)
@click.option(
    "-re",
    "--rescoring_engine",
    help="Either ristretto (MS²Rescore built-in) or percolator (run downstream) (default: `ristretto`)",
    type=click.Choice(RESCORING_ENGINES),
    default="ristretto",
)
@click.option("--train_fdr", help="FDR threshold for ristretto's semi-supervised training (default: `0.01`)", type=float, default=0.01)
@click.option("-d", "--id_decoy_pattern", help="Regex decoy pattern (default: `DECOY_`)", default="^DECOY_")
@click.option(
    "--fdr_level",
    help="FDR level written as `q-value` metavalue when rescoring with ristretto (default: `psm_level_fdrs`)",
    type=click.Choice(sorted(FDR_LEVELS)),
    default="psm_level_fdrs",
)
@click.option(
    "--require_precomputed_features",
    is_flag=True,
    default=False,
    help="Fail if any feature generator has to run, i.e. the input must already carry all MS²Rescore features "
    "(used for the dataset-wide ristretto pass on merged per-group output).",
)
def main(**kwargs):
    config = parse_cli_arguments_to_config(**kwargs)
    logging.info("MS²Rescore config:")
    logging.info(config)
    rescore_idxml(
        kwargs["psm_file"],
        kwargs["output_path"],
        config,
        kwargs["rescoring_engine"],
        kwargs["fdr_level"],
        kwargs["require_precomputed_features"],
    )


if __name__ == "__main__":
    sys.exit(main())
