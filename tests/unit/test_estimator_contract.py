"""Estimator contract (spec §8.4): ``fit`` returns ``self``; no function mutates its inputs."""

import copy

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2 import (
    DADA,
    BimeraRemover,
    BinnedQualErrorModel,
    ErrorLearner,
    LoessErrorModel,
    NoQualErrorModel,
    PacBioErrorModel,
    SpeciesAssigner,
    TaxonomyClassifier,
)
from pydada2.datasets import err_balanced_f, err_balanced_r, example_path

from rparity import load, matrix

TRAIN = example_path("example_train_set.fa.gz")
SPECIES = example_path("example_species_assignment.fa.gz")


@pytest.fixture(scope="module")
def dereps():
    return (
        dada2.derep_fastq(example_path("sam1F.fastq.gz")),
        dada2.derep_fastq(example_path("sam1R.fastq.gz")),
    )


def _snapshot_derep(d):
    return (d.uniques.copy(), d.quals.copy(), d.map.copy())


def _assert_derep_unchanged(d, snap):
    pd.testing.assert_series_equal(d.uniques, snap[0])
    np.testing.assert_array_equal(d.quals, snap[1])
    np.testing.assert_array_equal(d.map, snap[2])


def _trans():
    return matrix(load("error_models")["trans"])


@pytest.mark.parametrize(
    "est",
    [
        LoessErrorModel(),
        NoQualErrorModel(),
        PacBioErrorModel(),
        BinnedQualErrorModel(load("error_models")["binned_q"]),
    ],
    ids=lambda e: type(e).__name__,
)
def test_error_model_fit_returns_self_and_keeps_input(est):
    trans = _trans()
    if isinstance(est, PacBioErrorModel):
        trans = np.hstack([trans, np.zeros((16, 94 - trans.shape[1]))])
        trans[:, 93] = trans[:, 40]
    before = trans.copy()
    assert est.fit(trans) is est
    np.testing.assert_array_equal(trans, before)


def test_dada_fit_returns_self_and_keeps_inputs(dereps):
    drp = dereps[0]
    snap = _snapshot_derep(drp)
    err = err_balanced_f()[:, :30].copy()  # narrower than the quals: DADA extends a copy
    err_before = err.copy()
    priors = [drp.sequences[2]]
    est = DADA(err=err, priors=priors, params={"band_size": 32})
    assert est.fit(drp) is est
    _assert_derep_unchanged(drp, snap)
    np.testing.assert_array_equal(err, err_before)
    assert priors == [drp.sequences[2]]
    assert est.get_params()["params"] == {"band_size": 32}


def test_error_learner_fit_returns_self(dereps):
    drp = dereps[0]
    snap = _snapshot_derep(drp)
    est = ErrorLearner(n_bases=1e5)
    assert est.fit(drp) is est
    _assert_derep_unchanged(drp, snap)


def test_taxonomy_and_species_fit_return_self():
    clf = TaxonomyClassifier(random_state=1)
    assert clf.fit(TRAIN) is clf
    seqs = load("taxonomy")["queries"][:3]
    before = list(seqs)
    clf.predict(seqs)
    assert seqs == before
    sp = SpeciesAssigner()
    assert sp.fit(SPECIES) is sp
    sp.predict(seqs)
    assert seqs == before


def _seqtab():
    st = load("chimeras")["st_chim"]
    return pd.DataFrame(matrix(st, dtype=np.int64), index=st["rownames"], columns=st["colnames"])


@pytest.mark.parametrize("method", ["consensus", "pooled", "per-sample"])
def test_bimera_remover_fit_returns_self_and_keeps_table(method):
    st = _seqtab()
    before = st.copy()
    br = BimeraRemover(method=method)
    assert br.fit(st) is br
    br.transform(st)
    dada2.remove_bimera_denovo(st, method=method)
    dada2.is_bimera_denovo_table(st)
    dada2.is_bimera_denovo(st.sum(axis=0))
    pd.testing.assert_frame_equal(st, before)


def test_table_functions_keep_inputs():
    st = _seqtab()
    before = st.copy()
    dada2.collapse_no_mismatch(st)
    other = st.iloc[:, :3].rename(index=lambda s: s + "_b")
    other_before = other.copy()
    dada2.merge_sequence_tables(st, other)
    dada2.is_shift_denovo(st.sum(axis=0))
    pd.testing.assert_frame_equal(st, before)
    pd.testing.assert_frame_equal(other, other_before)
    samples = {"a": {"ACGT": 3, "GGCC": 1}, "b": {"ACGT": 1}}
    samples_before = copy.deepcopy(samples)
    dada2.make_sequence_table(samples)
    assert samples == samples_before


def test_merge_pairs_keeps_inputs(dereps):
    drp_f, drp_r = dereps
    dd_f = DADA(err=err_balanced_f()).fit(drp_f)[0]
    dd_r = DADA(err=err_balanced_r()).fit(drp_r)[0]
    snaps = [_snapshot_derep(drp_f), _snapshot_derep(drp_r)]
    maps = [dd_f.map.copy(), dd_r.map.copy()]
    clusters = [dd_f.clusters.copy(), dd_r.clusters.copy()]
    dada2.merge_pairs(dd_f, drp_f, dd_r, drp_r, return_rejects=True, propagate_col=["n0"])
    _assert_derep_unchanged(drp_f, snaps[0])
    _assert_derep_unchanged(drp_r, snaps[1])
    np.testing.assert_array_equal(dd_f.map, maps[0])
    np.testing.assert_array_equal(dd_r.map, maps[1])
    pd.testing.assert_frame_equal(dd_f.clusters, clusters[0])
    pd.testing.assert_frame_equal(dd_r.clusters, clusters[1])


def test_inflate_err_and_get_errors_keep_inputs():
    err = err_balanced_f()
    before = err.copy()
    dada2.inflate_err(err, 3, inflate_self_transitions=True)
    dada2.get_errors(err)
    np.testing.assert_array_equal(err, before)
