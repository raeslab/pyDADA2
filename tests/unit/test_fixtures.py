import json

import numpy as np

from rparity import DATA, load, matrix, scalar


def test_manifest():
    manifest = json.loads((DATA / "MANIFEST.json").read_text())
    assert manifest["dada2_version"] == "1.36.0"
    assert manifest["upstream_commit"] == "b61fe2e"


def test_namespace_exports():
    exports = (DATA / "NAMESPACE_exports.txt").read_text().split()
    assert len(exports) == 38
    assert "dada" in exports


def test_bundled_err_npy_matches_json():
    js = load("bundled_err")
    for name in ("errBalancedF", "errBalancedR", "tperr1"):
        arr = np.load(DATA / "err" / f"{name}.npy")
        assert arr.shape == (16, 41)
        np.testing.assert_array_equal(arr, matrix(js[name]))


def test_scalar():
    assert scalar([3]) == 3
    assert scalar(3) == 3
