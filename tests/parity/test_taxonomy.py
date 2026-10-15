"""Phase 8: taxonomy vs R (spec §8.3: exact top hits; statistical bootstrap agreement)."""

import numpy as np
import pandas as pd
import pytest

from pydada2.datasets import example_path
from pydada2.taxonomy import TaxonomyClassifier, add_species, assign_species

from rparity import load

pytestmark = pytest.mark.parity

TRAIN = example_path("example_train_set.fa.gz")
SPECIES = example_path("example_species_assignment.fa.gz")


@pytest.fixture(scope="module")
def tx():
    return load("taxonomy")


def _r_tax(obj):
    # Build the object array directly: pandas 3 would turn None into NaN in string columns.
    cols = sorted(obj, key=lambda c: int(c[1:]))
    return np.array([obj[c] for c in cols], dtype=object).T


def _py(df):
    return df.to_numpy(dtype=object)


def _tied(tx, key, i):
    return len(tx[key][i]) > 1


@pytest.mark.parametrize(
    "which,qkey,tkey",
    [("boot0", "queries", "tie_sets_queries"), ("boot0_full", "queries_full", "tie_sets_full")],
)
def test_top_hits_exact(tx, which, qkey, tkey):
    clf = TaxonomyClassifier(min_boot=0, random_state=100).fit(TRAIN)
    got = _py(clf.predict(tx[qkey]))
    exp = _r_tax(tx[which]["tax"])
    for i in range(len(exp)):
        if _tied(tx, tkey, i):
            assert got[i, -1] in tx[tkey][i]
        else:
            assert list(got[i]) == list(exp[i]), i


def test_try_rc_top_hits(tx):
    clf = TaxonomyClassifier(min_boot=0, try_rc=True, random_state=1).fit(TRAIN)
    got = _py(clf.predict(tx["tryrc"]["queries"]))
    np.testing.assert_array_equal(got, _r_tax(tx["tryrc"]["tax"]))


def _binomial_bound(p_a, p_b=None):
    """3 sigma of the difference of two independent 100-draw bootstrap estimates (+2 for ties),
    using the pooled estimate of the underlying proportion."""
    pooled = p_a if p_b is None else (np.asarray(p_a) + np.asarray(p_b)) / 2
    p = np.clip(np.asarray(pooled) / 100.0, 0.0, 1.0)
    return 3 * np.sqrt(2 * 100 * p * (1 - p)) + 2


@pytest.mark.parametrize("which,qkey", [("boot0", "queries"), ("boot0_full", "queries_full")])
def test_bootstrap_statistical(tx, which, qkey):
    # Bootstrap values differ from R only by RNG noise (DEVIATIONS.md #2): every entry within
    # 3 sigma of binomial noise, median |difference| <= 10 (spec §8.3).
    clf = TaxonomyClassifier(min_boot=0, random_state=100).fit(TRAIN)
    boot = clf.predict_bootstrap(tx[qkey]).to_numpy()
    exp = np.array(tx[which]["boot"])
    diff = np.abs(boot - exp)
    assert np.median(diff) <= 10
    assert (diff <= _binomial_bound(boot, exp)).all()


@pytest.mark.parametrize(
    "which,qkey,bkey",
    [("minboot50", "queries", "boot0"), ("minboot50_full", "queries_full", "boot0_full")],
)
def test_min_boot_50_agreement(tx, which, qkey, bkey):
    # Assignments must agree at every level whose confidence is clearly away from the threshold:
    # a disagreement is only allowed where R's bootstrap is within its binomial noise of 50.
    clf = TaxonomyClassifier(random_state=100).fit(TRAIN)
    got = _py(clf.predict(tx[qkey]))
    exp = _r_tax(tx[which])
    r_boot = np.array(tx[bkey]["boot"])
    for i in range(len(exp)):
        for j in range(exp.shape[1]):
            if got[i, j] == exp[i, j] or (got[i, j] is None and exp[i, j] is None):
                continue
            assert abs(r_boot[i, j] - 50) <= _binomial_bound(np.array(50.0)), (i, j)
            break  # deeper levels follow from this one


def test_random_state_and_threads(tx):
    q = tx["queries_full"]
    a = TaxonomyClassifier(random_state=7).fit(TRAIN)
    b = TaxonomyClassifier(random_state=7, n_jobs=4).fit(TRAIN)
    ta, ba = a.predict(q, return_bootstrap=True)
    tb, bb = b.predict(q, return_bootstrap=True)
    assert ta.equals(tb) and ba.equals(bb)


@pytest.mark.parametrize("case,kwargs", [
    ("default", {}), ("multiple", {"allow_multiple": True}), ("multiple3", {"allow_multiple": 3}),
])  # fmt: skip
def test_assign_species(tx, case, kwargs):
    got = _py(assign_species(tx["queries_species"], SPECIES, **kwargs))
    np.testing.assert_array_equal(got, _r_tax(tx["species"][case]))


def test_assign_species_try_rc_and_asv(tx):
    from pydada2.sequences import reverse_complement

    q = tx["queries_species"]
    mixed = q[:20] + reverse_complement(q[20:40])
    got = _py(assign_species(mixed, SPECIES, try_rc=True))
    np.testing.assert_array_equal(got, _r_tax(tx["species"]["tryrc"]))
    got = _py(assign_species(tx["queries"], SPECIES))
    np.testing.assert_array_equal(got, _r_tax(tx["species"]["asv"]))


@pytest.mark.parametrize(
    "case,kwargs", [("add_species", {}), ("add_species_multiple", {"allow_multiple": True})]
)
def test_add_species(tx, case, kwargs):
    r = _r_tax(tx[case])
    levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    taxtab = pd.DataFrame(r[:, :6], index=tx["queries_species"], columns=levels)
    got = add_species(taxtab, SPECIES, **kwargs)
    assert got["Species"].tolist() == list(r[:, 6])
