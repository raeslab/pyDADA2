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


def test_dada_uniques_validates_err_and_quals():
    seqs, ab, pri, err, quals = _dada_args()
    for bad in (np.nan, -0.1, 1.5, np.inf):
        e = err.copy()
        e[3, 7] = bad
        with pytest.raises(ValueError, match="Error matrix entries"):
            core.dada_uniques(seqs, ab, pri, e, quals)
    with pytest.raises(ValueError, match="Quality scores are required"):
        core.dada_uniques(seqs, ab, pri, err, None)
    # without quals and use_quals=False the transitions go to the first column (no crash)
    out = core.dada_uniques(seqs, ab, pri, np.ascontiguousarray(err[:, :1]), None, use_quals=False)
    assert out["trans"][:, 1:].sum() == 0 and out["trans"][:, 0].sum() > 0


def test_prior_born_clusters_have_birth_from():
    # omega_a = 0 forbids abundance births; with every sequence a prior and a lenient omega_p,
    # clusters are born from priors ("P"). birth_from was left uninitialised upstream.
    seqs, ab, _pri, err, quals = _dada_args()
    outs = [
        core.dada_uniques(seqs, ab, [True] * len(seqs), err, quals, omega_a=0.0, omega_p=0.5)
        for _ in range(3)
    ]
    bf = outs[0]["birth_from"]
    n = len(outs[0]["sequence"])
    assert n > 1
    assert bf[0] == -1 or bf[0] == 0
    assert ((bf[1:] >= 0) & (bf[1:] < n)).all()
    for o in outs[1:]:
        np.testing.assert_array_equal(o["birth_from"], bf)


def test_dada_uniques_raising_log_callback_propagates():
    seqs, ab, pri, err, quals = _dada_args()

    def boom(msg):
        raise RuntimeError("log failed")

    for _ in range(3):  # the B object is released on every abort (checked under ASan)
        with pytest.raises(RuntimeError, match="log failed"):
            core.dada_uniques(seqs, ab, pri, err, quals, verbose=True, log=boom)


def _samples():
    from pydada2 import derep_fastq
    from pydada2.datasets import err_balanced_f, example_path

    err = np.ascontiguousarray(err_balanced_f())
    out = []
    for name in ("sam1F.fastq.gz", "sam2F.fastq.gz"):
        drp = derep_fastq(example_path(name))
        for cut in (None, 400):  # a smaller sub-sample as a third/fourth sample
            seqs = drp.sequences[:cut]
            n = len(seqs)
            out.append(
                (seqs, drp.uniques.to_numpy()[:n].tolist(), [False] * n, err,
                 np.ascontiguousarray(drp.quals[:n]))
            )  # fmt: skip
    return out


def _batch(samples, **kw):
    cols = list(zip(*samples, strict=True))
    return core.dada_uniques_batch(*[list(c) for c in cols], **kw)


def _same(a, b):
    assert a["sequence"] == b["sequence"]
    for key in ("abundance", "map", "pval_raw", "trans", "birth_from", "birth_pval"):
        np.testing.assert_array_equal(a[key], b[key])


@pytest.mark.parametrize(("n_threads", "sample_threads"), [(1, 1), (4, 1), (1, 4), (4, 4)])
def test_dada_uniques_batch_equals_single_calls(n_threads, sample_threads):
    samples = _samples()
    ref = [core.dada_uniques(*s) for s in samples]
    outs = _batch(samples, n_threads=n_threads, sample_threads=sample_threads)
    assert len(outs) == len(ref)
    for a, b in zip(outs, ref, strict=True):
        _same(a, b)


@pytest.mark.parametrize("budget", [0, 20_000, 1 << 30])
def test_alignment_cache_gives_identical_results(budget):
    samples = _samples()
    ref = [core.dada_uniques(*s) for s in samples]
    cache = core.AlignmentCache(len(samples), budget)
    for _ in range(3):  # repeated rounds reuse the cached comparisons
        outs = _batch(samples, n_threads=2, sample_threads=2, cache=cache)
        for a, b in zip(outs, ref, strict=True):
            _same(a, b)
    assert cache.bytes_used <= budget
    if budget > 1 << 20:
        assert cache.bytes_used > 0
    # different samples in the same slots (the slots are reset), other alignment parameters
    outs = _batch(samples[::-1], cache=cache)
    for a, b in zip(outs, ref[::-1], strict=True):
        _same(a, b)
    ref_band = [core.dada_uniques(*s, band_size=8) for s in samples]
    for a, b in zip(_batch(samples, band_size=8, cache=cache), ref_band, strict=True):
        _same(a, b)


def test_alignment_cache_with_changing_error_rates():
    # the cache stores alignments only: rounds with other error matrices give fresh results
    samples = _samples()[:2]
    cache = core.AlignmentCache(2, 1 << 30)
    _batch(samples, cache=cache)
    err2 = np.ascontiguousarray(np.clip(samples[0][3] * 1.7, 0, 1))
    mod = [(s[0], s[1], s[2], err2, s[4]) for s in samples]
    ref = [core.dada_uniques(*s) for s in mod]
    for a, b in zip(_batch(mod, cache=cache, n_threads=3), ref, strict=True):
        _same(a, b)


def test_dada_uniques_batch_error_stops_batch():
    samples = _samples()
    bad = list(samples[1])
    bad[0] = ["ACGTN" * 50, *bad[0][1:]]
    with pytest.raises(ValueError, match="A/C/G/T"):
        _batch([samples[0], tuple(bad), samples[2]], sample_threads=3)
    with pytest.raises(ValueError, match="one entry per sample"):
        core.dada_uniques_batch([samples[0][0]], [], [], [], [])
    assert core.dada_uniques_batch([], [], [], [], []) == []


@pytest.mark.skipif(not hasattr(signal, "pthread_kill"), reason="needs pthread_kill")
def test_dada_uniques_batch_interrupt():
    samples = _samples() * 12
    main = threading.main_thread().ident
    timer = threading.Timer(0.05, signal.pthread_kill, args=(main, signal.SIGINT))
    timer.start()
    try:
        with pytest.raises(KeyboardInterrupt):
            _batch(samples, n_threads=2, sample_threads=2)
    finally:
        timer.cancel()
    # the pool and the handlers still work afterwards
    _same(_batch(samples[:1], n_threads=2)[0], core.dada_uniques(*samples[0]))


def _match_ref_reference(seqs, ref, w, non_overlapping):
    # upstream C_matchRef semantics: words of the circularised reference, scan with skip
    ref2 = ref + ref[:w]
    words = {ref2[i : i + w] for i in range(len(ref))}
    out = []
    for s in seqs:
        count, j = 0, 0
        while len(s) >= w and j <= len(s) - w:
            if s[j : j + w] in words:
                count += 1
                if non_overlapping:
                    j += w
            j += 1
        out.append(count)
    return out


def test_match_ref_matches_upstream_semantics():
    from pydada2.datasets import phix_genome

    rng = np.random.default_rng(11)
    g = phix_genome()
    seqs = []
    for _ in range(600):
        n = int(rng.choice([3, 15, 16, 17, 60, 240]))
        if rng.random() < 0.6:
            st = int(rng.integers(0, len(g) - n))
            s = list(g[st : st + n])
        else:
            s = list(rng.choice(list("ACGT"), n))
        for _ in range(int(rng.choice([0, 0, 1, 4]))):
            if s:
                s[int(rng.integers(len(s)))] = str(rng.choice(list("ACGTNa-")))
        seqs.append("".join(s))
    for ref in (g[:900], "ACGTACGTAC", "ACGTNACGTACGTTT"):
        for w in (1, 5, 16, 32, 33):
            for nono in (True, False):
                exp = _match_ref_reference(seqs, ref, w, nono)
                for nt in (1, 4):
                    assert list(core.match_ref(seqs, ref, w, nono, nt)) == exp
