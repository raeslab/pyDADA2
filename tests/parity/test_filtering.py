"""Phase 3: filter_and_trim / remove_primers parity with R (output records via canonical md5)."""

import hashlib
import math
from pathlib import Path

import dnaio
import pytest

import pydada2 as dada2

from rparity import EXTDATA, INPUTS, load, scalar

pytestmark = pytest.mark.parity

R2PY = {
    "truncLen": "trunc_len", "maxN": "max_n", "maxEE": "max_ee", "truncQ": "trunc_q",
    "rm.phix": "rm_phix", "trimLeft": "trim_left", "trimRight": "trim_right", "minLen": "min_len",
    "maxLen": "max_len", "minQ": "min_q", "rm.lowcomplex": "rm_lowcomplex", "n": "n",
    "orient.fwd": "orient_fwd", "matchIDs": "match_ids",
}  # fmt: skip


def _digest(path: Path) -> dict:
    if not path.exists():
        return {"exists": False}
    ids, widths, txt = [], [], []
    with dnaio.open(path) as fh:
        for r in fh:
            ids.append(r.name)
            widths.append(len(r.sequence))
            txt.append(f"{r.name}\n{r.sequence}\n{r.qualities}\n")
    return {
        "exists": True,
        "n": len(ids),
        "ids": ids,
        "widths": widths,
        "md5": hashlib.md5("".join(txt).encode()).hexdigest(),
    }


def _input(name: str) -> Path:
    return EXTDATA / name if (EXTDATA / name).exists() else INPUTS / name


def _params(rparams: dict, paired: bool) -> dict:
    out = {}
    for k, v in rparams.items():
        if isinstance(v, list):
            v = [math.inf if x == "Inf" else x for x in v]
            v = v[0] if (len(v) == 1 or not paired) else tuple(v)
        out[R2PY[k]] = v
    if "n" in out:
        out["n"] = int(out["n"])
    return out


def _check_output(path: Path, exp: dict) -> None:
    got = _digest(path)
    assert got["exists"] == scalar(exp["exists"]), path
    if got["exists"]:
        assert got["n"] == scalar(exp["n"])
        assert got["ids"] == exp["ids"]
        assert got["widths"] == exp["widths"]
        assert got["md5"] == scalar(exp["md5"])


CASES = sorted(load("filter_and_trim").keys())


@pytest.mark.parametrize("case", CASES)
def test_filter_and_trim_parity(case, tmp_path):
    fx = load("filter_and_trim")[case]
    paired = fx.get("inputs_rev") is not None
    fwd = [_input(f) for f in fx["inputs"]]
    filt = [tmp_path / f"F_{i}_{p.name}" for i, p in enumerate(fwd)]
    kwargs = _params(fx["params"] if isinstance(fx["params"], dict) else {}, paired)
    rev = filt_rev = None
    if paired:
        rev = [_input(f) for f in fx["inputs_rev"]]
        filt_rev = [tmp_path / f"R_{i}_{p.name}" for i, p in enumerate(rev)]
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = dada2.filter_and_trim(fwd, filt, rev, filt_rev, compress=True, **kwargs)
    assert df["reads_in"].tolist() == fx["reads_in"]
    assert df["reads_out"].tolist() == fx["reads_out"]
    for p, exp in zip(filt, fx["outputs"], strict=True):
        _check_output(p, exp)
    if paired:
        for p, exp in zip(filt_rev, fx["outputs_rev"], strict=True):
            _check_output(p, exp)


def test_filter_n_jobs_identical(tmp_path):
    fwd = [EXTDATA / "sam1F.fastq.gz", EXTDATA / "sam2F.fastq.gz"]
    a = [tmp_path / "a1.fq", tmp_path / "a2.fq"]
    b = [tmp_path / "b1.fq", tmp_path / "b2.fq"]
    d1 = dada2.filter_and_trim(fwd, a, trunc_len=240, max_ee=2, compress=False)
    d2 = dada2.filter_and_trim(fwd, b, trunc_len=240, max_ee=2, compress=False, n_jobs=4)
    assert d1.equals(d2)
    for x, y in zip(a, b, strict=True):
        assert x.read_bytes() == y.read_bytes()


PRIMER_CASES = ["default", "indels", "fwd_only", "no_trim", "mm0"]
PRIMER_ARGS = {
    "primer.fwd": "primer_fwd", "primer.rev": "primer_rev", "orient": "orient",
    "allow.indels": "allow_indels", "trim.fwd": "trim_fwd", "trim.rev": "trim_rev",
    "max.mismatch": "max_mismatch",
}  # fmt: skip


@pytest.mark.parametrize("case", PRIMER_CASES)
def test_remove_primers_parity(case, tmp_path):
    fx = load("remove_primers")[case]
    kwargs = {PRIMER_ARGS[k]: scalar(v) for k, v in fx["params"].items()}
    fout = tmp_path / "out.fastq.gz"
    df = dada2.remove_primers(EXTDATA / "samPBprimers.fastq.gz", fout, **kwargs)
    assert df["reads_in"].tolist() == fx["reads_in"]
    assert df["reads_out"].tolist() == fx["reads_out"]
    _check_output(fout, fx["output"])
