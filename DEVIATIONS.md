# DEVIATIONS.md — known differences between pydada2 and R dada2 1.36

Every known difference from the R reference is listed here (spec §8.3). Each entry has a
root cause, the outputs it affects, its measured magnitude, and the test that pins that
magnitude.

| # | Area | Difference | Cause | Magnitude (example data) | Pinned by |
|---|---|---|---|---|---|
| 1 | ppois (abundance p-values) | Boost.Math `gamma_p` instead of R's nmath `ppois` | no GPL code (spec §5.4) | rel ≤ 1e-12 on R grid; subnormal outputs (< 1e-305) differ in the last few significant digits | `tests/parity/test_core_lowlevel.py::test_ppois_upper` |
| 3 | `quality_profile` | uses the first `n` reads | R's `qa()` samples `n` reads at random from larger files | identical on files with ≤ n reads (all fixtures) | `test_quality_profile_parity` |
| 4 | loess error model / self-consistent DADA | `skmisc.loess` agrees with R's `loess` to a few ULP (~1e-15 rel. in the fit, ~1e-14 in error rates after `10^x`) | different builds/versions of the netlib dloess code | error matrices rel ≤ 1e-12; in `self_consistent` runs, near-ties between *underflowed* (0 vs subnormal) abundance p-values can change the **order** in which clusters are created (cluster indices in `clusters`/`map`/`birth_subs`); ASVs, abundances and cluster membership identical on the example data | `test_dada.py::test_dada_parity[self_consist]` (order-insensitive), `test_error_models.py::test_loess` |
| 2 | taxonomy RNG | bootstrap uniforms from `std::mt19937_64(random_state)` instead of R's RNG; deterministic tie-breaking | spec §5.4 | top hits identical; bootstrap \|Δ\| median 0, max 15 (within binomial noise); at `min_boot=50` 7/109 sequences differ at one level, all with R bootstrap in 42–55 | `tests/parity/test_taxonomy.py` |

## Upstream non-determinism (not a pydada2 deviation)

- **Taxonomy tie-breaking.** R dada2 breaks ties between equally likely genera with
  `std::random_device`, so `set.seed()` does not make `assignTaxonomy` reproducible.
  Measured R-vs-R on the fixture data: bootstrap |Δ| ≤ 2 (median 0); assignments at
  `minBoot=50` 99/100 identical. pydada2 is reproducible for a given `random_state`.

## Documented R quirks preserved on purpose

- `fastqFilter`/`filterAndTrim`: `minQ` only applies when `minQ > truncQ`, and reads are kept
  when `min(qual) > minQ` (the docs say reads *below* minQ are discarded).

## Behaviour changes at the C++ level (no effect on valid inputs)

- Unequal-length aligned strings in `eval_pair`/`pair_consensus` raise `ValueError` (R: warning +
  `NULL`); empty strings give zero counts (R: out-of-bounds read).
- Upstream memory leak in `C_is_bimera` (alignment not freed on `continue`) fixed.
- Taxonomy: sequences longer than 9999 nt no longer overflow fixed-size stack buffers; queries
  shorter than the k-mer size no longer blow up the uniform-buffer size.
- On non-x86 platforms SSE modes are emulated exactly (R on ARM forces `SSE=0`, which can change
  the k-mer/gapless screens); pydada2 therefore matches R-on-x86 everywhere.
- Single-end `orient_fwd`: reads shorter than the orientation sequence are treated as not
  matching (R: `narrow()` error aborts the whole file).
- `rm_lowcomplex`: reads without any A/C/G/T k-mer (complexity NA in R) are removed (R would
  index with NA).
- `remove_primers` with an empty primer raises `ValueError`.

