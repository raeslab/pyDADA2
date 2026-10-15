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


def _dada_args():
    from pydada2 import derep_fastq
    from pydada2.datasets import err_balanced_f, example_path

    drp = derep_fastq(example_path("sam1F.fastq.gz"))
    seqs = drp.sequences
    return (
        seqs,
        drp.uniques.to_numpy().tolist(),
        [False] * len(seqs),
        np.ascontiguousarray(err_balanced_f()),
        np.ascontiguousarray(drp.quals),
    )


def test_concurrent_calls_with_logging_do_not_interfere():
    # Handlers are thread-local: concurrent verbose calls from several Python threads each log
    # through their own callback and give the same results as a serial call.
    from concurrent.futures import ThreadPoolExecutor

    args = _dada_args()
    ref = core.dada_uniques(*args)
    logs = {k: [] for k in range(4)}

    def run(k):
        return core.dada_uniques(*args, n_threads=2, verbose=True, log=logs[k].append)

    with ThreadPoolExecutor(4) as ex:
        outs = list(ex.map(run, range(4)))
    for k, out in enumerate(outs):
        assert out["sequence"] == ref["sequence"]
        np.testing.assert_array_equal(out["map"], ref["map"])
        assert any("ALIGN" in m for m in logs[k])
    assert all(logs[k] == logs[0] for k in logs)


def test_concurrent_assign_taxonomy():
    from concurrent.futures import ThreadPoolExecutor

    model, genera = _toy_model()
    seqs = [g[20:280] for g in genera]
    rcs = [_rc(s) for s in seqs]
    ref = core.assign_taxonomy(model, seqs, rcs, True, 7, 1)
    with ThreadPoolExecutor(8) as ex:
        outs = list(ex.map(lambda _: core.assign_taxonomy(model, seqs, rcs, True, 7, 2), range(32)))
    for out in outs:
        for key in ("tax", "boot", "boot_tax"):
            np.testing.assert_array_equal(out[key], ref[key])


def test_scoped_handlers_restore_previous():
    # A nested call that installs no log handler must not remove the outer call's handler:
    # after the inner (handler-less) call returns, the outer thread logs again.
    args = _dada_args()
    msgs = []
    core.dada_uniques(*args, verbose=True, log=msgs.append)
    n = len(msgs)
    core.dada_uniques(*args)
    core.dada_uniques(*args, verbose=True, log=msgs.append)
    assert len(msgs) == 2 * n


@pytest.mark.parametrize("n_threads", [2, 8])
def test_parallel_for_nested_regions(n_threads):
    vals, chunks = core._parallel_selftest(200, 3, n_threads, -1, 17)
    np.testing.assert_array_equal(vals, np.arange(200, dtype=np.int64) ** 2)
    assert chunks == -(-200 // 3)


def test_parallel_for_concurrent_callers_share_the_pool():
    from concurrent.futures import ThreadPoolExecutor

    def run(k):
        return core._parallel_selftest(5000 + k, 7, 4, -1, 3 if k % 2 else 0)

    with ThreadPoolExecutor(8) as ex:
        outs = list(ex.map(run, range(24)))
    for k, (vals, _) in enumerate(outs):
        np.testing.assert_array_equal(vals, np.arange(5000 + k, dtype=np.int64) ** 2)


def test_parallel_for_huge_thread_request_is_capped():
    vals, chunks = core._parallel_selftest(100, 1, 100_000)
    np.testing.assert_array_equal(vals, np.arange(100, dtype=np.int64) ** 2)
    assert chunks == 100


def test_nwvec_scratch_reuse_and_long_sequences():
    # The vectorised aligner reuses per-thread buffers: alignments must not depend on what was
    # aligned before (sizes growing and shrinking), and sequences longer than the old fixed
    # buffers (SEQLEN = 9999) must work.
    rng = np.random.default_rng(5)

    def pair(n):
        a = rng.choice(list("ACGT"), n)
        b = a.copy()
        idx = rng.choice(n, max(1, n // 50), replace=False)
        b[idx] = rng.choice(list("ACGT"), idx.size)
        return "".join(a), "".join(np.delete(b, idx[:2]))

    pairs = [pair(n) for n in (250, 12000, 40, 251, 3000, 250)]
    s1 = [p[0] for p in pairs]
    s2 = [p[1] for p in pairs]
    batch = core.nwvec(s1, s2, 5, -4, -8, 16, True)
    for (a, b), (al0, al1) in zip(pairs, batch, strict=True):
        assert al0.replace("-", "") == a and al1.replace("-", "") == b
        assert core.nwvec([a], [b], 5, -4, -8, 16, True)[0] == (al0, al1)
