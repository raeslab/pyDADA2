import logging

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2 import NotFittedError, SpeciesAssigner, TaxonomyClassifier
from pydada2.datasets import example_path
from pydada2.io import read_fasta
from pydada2.taxonomy import _map_hits, _match_genera, _r_strsplit, _seed

TRAIN = example_path("example_train_set.fa.gz")


def test_r_strsplit():
    assert _r_strsplit("a;b;", ";") == ["a", "b"]
    assert _r_strsplit("a;b;;", ";") == ["a", "b", ""]
    assert _r_strsplit("a", ";") == ["a"]


def test_seed():
    assert _seed(5) == 5
    assert isinstance(_seed(None), int)
    g = np.random.default_rng(1)
    assert _seed(g) == _seed(np.random.default_rng(1))


def test_map_hits_semantics():
    names = ["A", "B", None, "Escherichia coli", "Shigella"]
    assert _map_hits([0], names, 1) == "A"
    assert _map_hits([0, 1], names, 1) is None
    assert _map_hits([1, 0], names, 2) == "A/B"
    assert _map_hits([0, 2], names, 1) is None  # NA counts as a distinct value
    assert _map_hits([0, 2], names, 2) == "A"  # ... but is dropped when pasting
    assert _map_hits([2], names, 1) is None
    assert _map_hits([], names, 1) is None
    assert _map_hits([3, 4], names, 1) == "Escherichia/Shigella"


def test_match_genera():
    assert _match_genera("Bacillus", "Bacillus")
    assert _match_genera("Bacillus_1", "Bacillus")
    assert _match_genera("Escherichia/Shigella", "Shigella")
    assert not _match_genera("Bacillus", "Bac")
    assert not _match_genera(None, "Bacillus")
    assert not _match_genera("", "Bacillus")
    assert not _match_genera(float("nan"), "Bacillus")
    assert not _match_genera("Foo(bar)_1", "Foo(bar)")  # valid regex (groups), as in R
    assert _match_genera("Foo[_1", "Foo[")  # invalid regex: literal fallback (R would error)


def test_classifier_errors_and_levels(tmp_path):
    with pytest.raises(NotFittedError):
        TaxonomyClassifier().predict(["ACGT" * 20])
    bad = [("Bacteria Firmicutes Bacillus", "ACGT" * 20)]
    with pytest.raises(ValueError, match="assign_species"):
        TaxonomyClassifier().fit(bad)
    with pytest.raises(ValueError, match="Incorrect reference"):
        TaxonomyClassifier().fit([("Bacteria", "ACGT" * 20)])
    with pytest.warns(UserWarning, match="too short"), pytest.raises(ValueError, match="No usable"):
        TaxonomyClassifier().fit([("A;B;", "ACGT")])
    refs = read_fasta(TRAIN)
    clf = TaxonomyClassifier(tax_levels=("K", "P"), random_state=0).fit(refs)
    assert clf._columns() == ["K", "P", "Level3", "Level4", "Level5", "Level6"]
    with pytest.warns(UserWarning, match="shorter than 50"):
        out = clf.predict(["ACGT" * 5, refs[0][1]])
    assert out.iloc[0].isna().all()
    with pytest.raises(ValueError, match="No sequences"):
        clf.predict([])


def test_unite_parsing(caplog):
    seq = read_fasta(TRAIN)[0][1]
    recs = [
        (
            f"SH{i}|FU|re{'p' if i % 2 else 'f'}s|x|"
            f"k__Fungi;p__Asco;c__unidentified;o__O;f__F;g__G{i};s__G{i}_sp",
            seq,
        )
        for i in range(12)
    ]
    with caplog.at_level(logging.INFO, logger="pydada2"):
        clf = TaxonomyClassifier(random_state=0, verbose=1).fit(recs)
    assert any("UNITE" in m for m in caplog.messages)
    assert clf.n_levels_ == 7
    assert "_DADA2_UNSPECIFIED" in clf.genera_[0]


def test_species_assigner_errors():
    with pytest.raises(NotFittedError):
        SpeciesAssigner().predict(["ACGT"])
    with pytest.raises(ValueError, match="assign_taxonomy"):
        SpeciesAssigner().fit([("Bacteria;Firmicutes;Bacilli;Bacillales;", "ACGT")])
    with pytest.raises(ValueError, match="Incorrect reference"):
        SpeciesAssigner().fit([("onlyid", "ACGT")])
    with pytest.raises(ValueError, match="Empty"):
        SpeciesAssigner().fit([])
    sa = SpeciesAssigner(verbose=1).fit([("x Genus species", "AACCGGTTAACCGGTT")])
    with pytest.raises(ValueError, match="Non-ACGT"):
        sa.predict(["ACGN"])
    assert sa.predict(["CCGGTT"]).loc["CCGGTT", "Species"] == "species"


def test_add_species_without_genus_column():
    tab = pd.DataFrame({"Level1": ["Genus"]}, index=["CCGGTT"])
    out = dada2.add_species(tab, [("x Genus species", "AACCGGTTAACCGGTT")], verbose=1)
    assert out["Species"].tolist() == ["species"]
