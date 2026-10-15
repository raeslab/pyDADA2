"""Every R NAMESPACE export maps to a Python object (or is explicitly not ported)."""

import importlib

from rparity import DATA

R_TO_PY = {
    "addSpecies": "add_species",
    "assignSpecies": "assign_species",
    "assignTaxonomy": "assign_taxonomy",
    "collapseNoMismatch": "collapse_no_mismatch",
    "dada": "DADA",
    "derepFastq": "derep_fastq",
    "fastqFilter": "fastq_filter",
    "fastqPairedFilter": "fastq_paired_filter",
    "filterAndTrim": "filter_and_trim",
    "getDadaOpt": "DadaParams",
    "getErrors": "get_errors",
    "getSequences": "get_sequences",
    "getUniques": "get_uniques",
    "inflateErr": "inflate_err",
    "isBimera": "is_bimera",
    "isBimeraDenovo": "is_bimera_denovo",
    "isBimeraDenovoTable": "is_bimera_denovo_table",
    "isPhiX": "is_phix",
    "isShiftDenovo": "is_shift_denovo",
    "learnErrors": "learn_errors",
    "loessErrfun": "LoessErrorModel",
    "makeBinnedQualErrfun": "BinnedQualErrorModel",
    "makeSequenceTable": "make_sequence_table",
    "mergePairs": "merge_pairs",
    "mergeSequenceTables": "merge_sequence_tables",
    "noqualErrfun": "NoQualErrorModel",
    "nwalign": "align.nwalign",
    "nwhamming": "align.nwhamming",
    "PacBioErrfun": "PacBioErrorModel",
    "plotComplexity": "plotting.plot_complexity",
    "plotErrors": "plotting.plot_errors",
    "plotQualityProfile": "plotting.plot_quality_profile",
    "rc": "reverse_complement",
    "removeBimeraDenovo": "remove_bimera_denovo",
    "removePrimers": "remove_primers",
    "seqComplexity": "seq_complexity",
    "setDadaOpt": "DadaParams",
    "uniquesToFasta": "uniques_to_fasta",
}

NOT_PORTED: dict[str, str] = {}


def _resolve(path):
    mod = importlib.import_module("pydada2")
    *mods, attr = path.split(".")
    for m in mods:
        mod = importlib.import_module(f"{mod.__name__}.{m}")
    return getattr(mod, attr)


def test_every_export_is_covered():
    exports = (DATA / "NAMESPACE_exports.txt").read_text().split()
    missing = [e for e in exports if e not in R_TO_PY and e not in NOT_PORTED]
    assert not missing, f"R exports without a Python mapping: {missing}"


def test_mapped_objects_exist():
    for r_name, py in R_TO_PY.items():
        assert callable(_resolve(py)), r_name
