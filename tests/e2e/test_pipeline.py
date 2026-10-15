"""End-to-end MiSeq-style pipeline on the bundled example data vs the R pipeline (e2e.json)."""

import warnings

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2.datasets import err_balanced_f, err_balanced_r, example_path

from rparity import load, matrix

pytestmark = [pytest.mark.slow, pytest.mark.parity]


def test_vignette_pipeline(tmp_path):
    e2e = load("e2e")
    fwd = [example_path("sam1F.fastq.gz"), example_path("sam2F.fastq.gz")]
    rev = [example_path("sam1R.fastq.gz"), example_path("sam2R.fastq.gz")]
    filt_f = [tmp_path / "sam1F.fq.gz", tmp_path / "sam2F.fq.gz"]
    filt_r = [tmp_path / "sam1R.fq.gz", tmp_path / "sam2R.fq.gz"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = dada2.filter_and_trim(fwd, filt_f, rev, filt_r, trunc_len=(240, 160), max_n=0,
                                    max_ee=(2, 2), trunc_q=2, rm_phix=True, n_jobs=2)  # fmt: skip
    names = ["sam1", "sam2"]
    drp_f = dict(zip(names, dada2.derep_fastq(filt_f), strict=True))
    drp_r = dict(zip(names, dada2.derep_fastq(filt_r), strict=True))
    dd_f = dada2.DADA(err=err_balanced_f(), n_jobs=2).fit(drp_f)
    dd_r = dada2.DADA(err=err_balanced_r(), n_jobs=2).fit(drp_r)
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
    np.testing.assert_array_equal(nochim.to_numpy(), matrix(exp, dtype=np.int64))

    taxa = dada2.assign_taxonomy(nochim.columns, example_path("example_train_set.fa.gz"),
                                 random_state=100, n_jobs=2)  # fmt: skip
    obj = e2e["taxa"]
    exp_tax = np.array([obj[c] for c in sorted(obj, key=lambda c: int(c[1:]))], dtype=object).T
    # Kingdom and Phylum have bootstrap 100 here: identical to R. Deeper levels may differ where
    # the bootstrap is near min_boot (RNG differs; see tests/parity/test_taxonomy.py).
    np.testing.assert_array_equal(taxa.to_numpy(dtype=object)[:, :2], exp_tax[:, :2])
