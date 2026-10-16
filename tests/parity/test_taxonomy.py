"""Taxonomy vs R (exact top hits; statistical bootstrap agreement)."""

import numpy as np
import pandas as pd
import pytest

from pydada2.datasets import example_path
from pydada2.taxonomy import TaxonomyClassifier, add_species, assign_species

from rparity import (
    assert_bootstrap_agreement,
    load,
    r_taxonomy,
    taxonomy_disagreements,
)

pytestmark = pytest.mark.parity

TRAIN = example_path("example_train_set.fa.gz")
SPECIES = example_path("example_species_assignment.fa.gz")


@pytest.fixture(scope="module")
def tx():
    return load("taxonomy")


_r_tax = r_taxonomy


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


@pytest.mark.parametrize("which,qkey", [("boot0", "queries"), ("boot0_full", "queries_full")])
def test_bootstrap_statistical(tx, which, qkey):
    # Bootstrap values differ from R only by RNG noise (deviation #2): every entry within
    # 3 sigma of binomial noise, median |difference| <= 10 (docs/porting-notes.md, §2).
    clf = TaxonomyClassifier(min_boot=0, random_state=100).fit(TRAIN)
    boot = clf.predict_bootstrap(tx[qkey]).to_numpy()
    assert_bootstrap_agreement(boot, np.array(tx[which]["boot"]))


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
    taxonomy_disagreements(got, exp, np.array(tx[bkey]["boot"]))


def test_deviation_2_magnitude(tx):
    # Pins deviation #2 (docs/porting-notes.md) tighter than the statistical criterion, so it
    # cannot grow silently. Over the 109 queries: bootstrap |delta| median 0, max 15; 8 sequences
    # differ at min_boot=50 (7-9 depending on the R run: R's own tie-breaking is not
    # reproducible), all with R bootstrap 42-55. The bounds leave room for that R-side
    # variation only.
    deltas, disagreements = [], []
    for which, qkey, bkey in [
        ("minboot50", "queries", "boot0"),
        ("minboot50_full", "queries_full", "boot0_full"),
    ]:
        clf = TaxonomyClassifier(random_state=100).fit(TRAIN)
        tax, boot = clf.predict(tx[qkey], return_bootstrap=True)
        r_boot = np.array(tx[bkey]["boot"])
        deltas.append(np.abs(boot.to_numpy() - r_boot).ravel())
        disagreements += taxonomy_disagreements(_py(tax), _r_tax(tx[which]), r_boot)
    d = np.concatenate(deltas)
    assert np.median(d) == 0
    assert d.max() <= 20
    assert len(disagreements) <= 12
    assert all(abs(b - 50) <= 10 for _, _, b in disagreements), disagreements


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
