"""End-to-end vignette pipeline on the bundled example data vs the R pipeline (e2e.json).

filter_and_trim -> learn_errors -> derep -> DADA -> merge_pairs -> make_sequence_table ->
consensus chimera removal -> assign_taxonomy, as in the DADA2 tutorial. The track table and
seqtab.nochim must be identical to R. Taxonomy is compared at every level with the spec §8.3
criterion against R's results for the same ASVs (taxonomy.json `queries`; R's taxonomy is not
reproducible run to run, so it is not stored in e2e.json).
"""

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2.datasets import example_path

from rparity import (
    assert_bootstrap_agreement,
    load,
    matrix,
    r_taxonomy,
    r_warnings,
    taxonomy_disagreements,
)

pytestmark = [pytest.mark.slow, pytest.mark.parity]


def test_vignette_pipeline(tmp_path):
    e2e = load("e2e")
    fwd = [example_path("sam1F.fastq.gz"), example_path("sam2F.fastq.gz")]
    rev = [example_path("sam1R.fastq.gz"), example_path("sam2R.fastq.gz")]
    filt_f = [tmp_path / "sam1F.fq.gz", tmp_path / "sam2F.fq.gz"]
    filt_r = [tmp_path / "sam1R.fq.gz", tmp_path / "sam2R.fq.gz"]
    with r_warnings([]):  # R's filterAndTrim warns about nothing here
        out = dada2.filter_and_trim(fwd, filt_f, rev, filt_r, trunc_len=(240, 160), max_n=0,
                                    max_ee=(2, 2), trunc_q=2, rm_phix=True, n_jobs=2)  # fmt: skip
    err_f = dada2.learn_errors(filt_f, n_jobs=2)
    err_r = dada2.learn_errors(filt_r, n_jobs=2)
    names = ["sam1", "sam2"]
    drp_f = dict(zip(names, dada2.derep_fastq(filt_f), strict=True))
    drp_r = dict(zip(names, dada2.derep_fastq(filt_r), strict=True))
    dd_f = dada2.DADA(err=err_f, n_jobs=2).fit(drp_f)
    dd_r = dada2.DADA(err=err_r, n_jobs=2).fit(drp_r)
    mergers = dada2.merge_pairs(dd_f, list(drp_f.values()), dd_r, list(drp_r.values()))
    seqtab = dada2.make_sequence_table(dict(zip(names, mergers, strict=True)))
    nochim = dada2.BimeraRemover(method="consensus", n_jobs=2).fit_transform(seqtab)

    track = pd.DataFrame(
        {
            "input": out["reads_in"].to_numpy(),
            "filtered": out["reads_out"].to_numpy(),
            "denoisedF": [r.denoised.sum() for r in dd_f.results_],
            "denoisedR": [r.denoised.sum() for r in dd_r.results_],
            "merged": [m["abundance"].sum() for m in mergers],
            "nonchim": nochim.sum(axis=1).to_numpy(),
        },
        index=names,
    )
    np.testing.assert_array_equal(track.to_numpy(), matrix(e2e["track"], dtype=np.int64))

    exp = e2e["seqtab_nochim"]
    assert nochim.columns.tolist() == exp["colnames"]
    assert nochim.index.tolist() == exp["rownames"]
    np.testing.assert_array_equal(nochim.to_numpy(), matrix(exp, dtype=np.int64))

    tax, boot = dada2.assign_taxonomy(
        nochim.columns, example_path("example_train_set.fa.gz"), output_bootstraps=True,
        random_state=100, n_jobs=2,
    )  # fmt: skip
    tx = load("taxonomy")
    rows = [tx["queries"].index(s) for s in nochim.columns]
    r_tax = r_taxonomy(tx["minboot50"])[rows]
    r_boot = np.array(tx["boot0"]["boot"])[rows]
    assert_bootstrap_agreement(boot.to_numpy(), r_boot)
    taxonomy_disagreements(tax.to_numpy(dtype=object), r_tax, r_boot)
