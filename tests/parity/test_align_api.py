"""Public alignment API (pydada2.align) vs R's nwalign / nwhamming / C_eval_pair / C_pair_consensus.

tests/parity/test_core_lowlevel.py checks the bindings (pydada2._core); this file checks the
public wrappers with R's argument names and defaults on the same align.json panel.
"""

import numpy as np
import pytest

from pydada2 import align

from rparity import load, scalar

pytestmark = pytest.mark.parity

AL = load("align")
PANEL = AL["panel"]
PAIRS = np.array(AL["pairs"]) - 1

# R nwalign() arguments of each variant -> pydada2.align.nwalign keyword arguments.
R2PY = {"band": "band", "endsfree": "endsfree", "homo_gap": "homo_gap", "vec": "vec",
        "match": "match", "mismatch": "mismatch", "gap": "gap"}  # fmt: skip


def _kwargs(variant):
    v = AL["variants"][variant]
    if not isinstance(v, dict):  # empty R list -> []: R defaults
        return {}
    return {R2PY[k]: scalar(x) for k, x in v.items()}


@pytest.mark.parametrize("variant", sorted(AL["variants"]))
def test_nwalign_public_api(variant):
    res = AL["results"][variant]
    kwargs = _kwargs(variant)
    for k, (i, j) in enumerate(PAIRS):
        al1, al2 = align.nwalign(PANEL[i], PANEL[j], **kwargs)
        assert (al1, al2) == (res["al1"][k], res["al2"][k]), (variant, k)
        assert tuple(align.eval_pair(al1, al2)) == (
            res["match"][k], res["mismatch"][k], res["indel"][k],
        )  # fmt: skip


@pytest.mark.parametrize("key,kwargs", [("nwhamming", {}), ("nwhamming_band16", {"band": 16})])
def test_nwhamming_public_api(key, kwargs):
    s1 = [PANEL[i] for i in PAIRS[:, 0]]
    s2 = [PANEL[j] for j in PAIRS[:, 1]]
    got = align.nwhamming(s1, s2, **kwargs)
    assert np.asarray(got).tolist() == AL[key]
    # scalar form
    assert align.nwhamming(s1[0], s2[0], **kwargs) == AL[key][0]


def test_pair_consensus_public_api():
    for r in AL["pair_consensus"]["results"]:
        a1, a2 = scalar(r["al1"]), scalar(r["al2"])
        fseq = AL["pair_consensus"]["fwd"][scalar(r["i"]) - 1]
        rseq = AL["pair_consensus"]["rev_rc"][scalar(r["j"]) - 1]
        assert align.nwalign(fseq, rseq, match=1, mismatch=-8, gap=-8, band=-1) == (a1, a2)
        assert list(align.eval_pair(a1, a2)) == r["eval"]
        assert align.pair_consensus(a1, a2) == scalar(r["p1"])
        assert align.pair_consensus(a1, a2, prefer=2) == scalar(r["p2"])
        assert align.pair_consensus(a1, a2, trim_overhang=True) == scalar(r["p1_trim"])
        assert align.pair_consensus(a1, a2, prefer=2, trim_overhang=True) == scalar(r["p2_trim"])
