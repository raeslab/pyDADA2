import signal
import threading

import numpy as np
import pytest

import pydada2._core as core


@pytest.mark.parametrize("n_threads", [1, 2, 4, 16])
@pytest.mark.parametrize("grain", [1, 7, 1000])
def test_parallel_for_covers_range(n_threads, grain):
    vals, chunks = core._parallel_selftest(1000, grain, n_threads)
    np.testing.assert_array_equal(vals, np.arange(1000, dtype=np.int64) ** 2)
    assert chunks == (1 if n_threads == 1 else -(-1000 // grain))


def test_parallel_for_empty():
    vals, chunks = core._parallel_selftest(0, 10, 4)
    assert vals.size == 0 and chunks == 0


@pytest.mark.parametrize("n_threads", [1, 4])
def test_parallel_for_propagates_exceptions(n_threads):
    with pytest.raises(RuntimeError, match="parallel failure at 517"):
        core._parallel_selftest(1000, 3, n_threads, 517)


def test_interrupt_check_without_signal():
    core._interrupt_selftest(10)


@pytest.mark.skipif(not hasattr(signal, "pthread_kill"), reason="needs pthread_kill")
def test_interrupt_raises_keyboard_interrupt():
    # Deliver SIGINT while the C++ loop runs; the next interrupt check must raise KeyboardInterrupt.
    main = threading.main_thread().ident
    timer = threading.Timer(0.05, signal.pthread_kill, args=(main, signal.SIGINT))
    timer.start()
    try:
        with pytest.raises(KeyboardInterrupt):
            core._interrupt_selftest(10**12)
    finally:
        timer.cancel()


def _toy_model():
    rng = np.random.default_rng(0)
    bases = np.array(list("ACGT"))
    genera = ["".join(rng.choice(bases, 300)) for _ in range(6)]
    refs, r2g = [], []
    for g, base in enumerate(genera):
        for _ in range(4):
            s = np.array(list(base))
            idx = rng.choice(300, 6, replace=False)
            s[idx] = rng.choice(bases, 6)
            refs.append("".join(s))
            r2g.append(g)
    genusmat = np.array([[0, g // 3, g] for g in range(6)], dtype=np.int32)
    return core.TaxonomyModel(refs, r2g, genusmat), genera


def _rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def test_taxonomy_smoke_and_determinism():
    model, genera = _toy_model()
    assert (model.ngenus, model.nlevel) == (6, 3)
    seqs = [g[20:280] for g in genera] + ["ACGT" * 10]  # last one too short (< 50 nt)
    rcs = [_rc(s) for s in seqs]
    r1 = core.assign_taxonomy(model, seqs, rcs, False, 42, 1)
    assert r1["tax"].tolist() == [0, 1, 2, 3, 4, 5, -1]
    assert (r1["boot"][:6, 2] >= 90).all()
    assert (r1["boot_tax"][6] == -1).all() and (r1["boot"][6] == 0).all()
    for n in (2, 4):
        rn = core.assign_taxonomy(model, seqs, rcs, False, 42, n)
        for key in ("tax", "boot", "boot_tax"):
            np.testing.assert_array_equal(rn[key], r1[key])
    # reverse-complemented queries are recovered with try_rc
    rr = core.assign_taxonomy(model, rcs[:6], seqs[:6], True, 42, 1)
    assert rr["tax"].tolist() == [0, 1, 2, 3, 4, 5]


def test_taxonomy_validation():
    with pytest.raises(ValueError, match="No reference"):
        core.TaxonomyModel([], [], np.zeros((1, 1), dtype=np.int32))
    with pytest.raises(ValueError, match="Invalid map"):
        core.TaxonomyModel(["ACGTACGTACGT"], [3], np.zeros((1, 1), dtype=np.int32))
    model, _ = _toy_model()
    with pytest.raises(ValueError, match="No seqs"):
        core.assign_taxonomy(model, [], [], False, 1, 1)
