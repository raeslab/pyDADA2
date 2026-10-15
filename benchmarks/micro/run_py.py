"""Micro-benchmarks of pydada2 on the bundled example data (spec §9.1).

Usage: uv run python benchmarks/micro/run_py.py [repeats] [outfile]

Stage names, inputs and the timing method mirror benchmarks/micro/run_r.R:

* every stage runs once as warm-up; stages faster than ``MIN_SAMPLE`` seconds are then timed in
  loops of k calls (k chosen from the warm-up), and the reported time is the median over
  ``repeats`` samples of (loop time / k);
* ``gc.collect()`` runs before every timed sample (outside the timed region), like R's gc();
* thread counts equal R's: filter_and_trim at 1 and 2 (R parallelises over the 2 files only),
  the multithreaded R stages at 1 and 4. Stages R cannot run multithreaded are additionally
  timed at 4 threads for pydada2 only (keys ending in "@4" without an R counterpart);
* inputs are prepared outside the timed regions.
"""

from __future__ import annotations

import gc
import json
import math
import os
import platform
import statistics
import sys
import tempfile
import time
import warnings
from collections.abc import Callable
from pathlib import Path

import pydada2 as dada2
from pydada2.datasets import err_balanced_f, err_balanced_r, example_path
from pydada2.io import read_fasta

MIN_SAMPLE = 0.2  # seconds
LOOPS: dict[str, int] = {}


def bench(key: str, fn: Callable[[], object], repeats: int) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gc.collect()
        t0 = time.perf_counter()
        fn()
        warm = time.perf_counter() - t0
        k = int(max(1, min(1000, math.ceil(MIN_SAMPLE / max(warm, 1e-6)))))
        times = []
        for _ in range(repeats):
            gc.collect()
            t0 = time.perf_counter()
            for _ in range(k):
                fn()
            times.append((time.perf_counter() - t0) / k)
    LOOPS[key] = k
    return statistics.median(times)


def chimera_table(st):  # type: ignore[no-untyped-def]
    """st plus chimeras of its 6 most abundant sequences (breakpoint 120), as in run_r.R."""
    import pandas as pd

    top = list(st.columns[:6])
    chims = [a[:120] + b[120:] for a in top for b in top if a != b]
    chims = [c for c in dict.fromkeys(chims) if c not in set(st.columns)]
    extra = pd.DataFrame({c: [5, 3] for c in chims}, index=st.index)
    return pd.concat([st, extra], axis=1), len(chims)


def main() -> None:
    repeats = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "benchmarks/micro/results_py.json")
    tmp = Path(tempfile.mkdtemp())
    fwd = [example_path("sam1F.fastq.gz"), example_path("sam2F.fastq.gz")]
    rev = [example_path("sam1R.fastq.gz"), example_path("sam2R.fastq.gz")]
    filt_f = [tmp / "F1.fq.gz", tmp / "F2.fq.gz"]
    filt_r = [tmp / "R1.fq.gz", tmp / "R2.fq.gz"]
    ef, er = err_balanced_f(), err_balanced_r()
    filt_kw = dict(trunc_len=(240, 160), max_n=0, max_ee=(2, 2), trunc_q=2, rm_phix=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dada2.filter_and_trim(fwd, filt_f, rev, filt_r, **filt_kw)
    drp_f = dada2.derep_fastq(filt_f)
    drp_r = dada2.derep_fastq(filt_r)
    dd_f = dada2.DADA(err=ef).fit(drp_f)
    dd_r = dada2.DADA(err=er).fit(drp_r)
    mergers = dada2.merge_pairs(dd_f, drp_f, dd_r, drp_r)
    st = dada2.make_sequence_table(dict(zip(["sam1", "sam2"], mergers, strict=True)))
    st_chim, n_chims = chimera_table(st)
    ten = [s for _, s in read_fasta(example_path("ten_16s.100.fa.gz"))]
    qfull = ten[::10]
    species_ref = example_path("example_species_assignment.fa.gz")
    qspec = [s[i : i + 250] for _, s in read_fasta(species_ref) for i in range(0, len(s) - 249, 50)]
    qspec = [s for s in qspec if set(s) <= set("ACGT")]
    f27, r1492rc = "AGRGTTYGATYMTGGCTCAG", dada2.reverse_complement("RGYTACCTTGTTACGACTT")
    train = example_path("example_train_set.fa.gz")
    sam1_seqs = dada2.get_sequences(example_path("sam1F.fastq.gz"))
    spec_hits = dada2.assign_species(qspec, species_ref)

    res: dict[str, object] = {
        "meta": {
            "python": platform.python_version(),
            "pydada2": dada2.__version__,
            "repeats": repeats,
            "cores": os.cpu_count(),
            "min_sample": MIN_SAMPLE,
            "n_chimeras": n_chims,
            "n_species_queries": len(qspec),
            "species_hits": int(spec_hits["Species"].notna().sum()),
        }
    }

    def run(key: str, fn: Callable[[], object]) -> None:
        res[key] = bench(key, fn, repeats)

    for nt in (1, 2):
        run(
            f"filter_and_trim@{nt}",
            lambda nt=nt: dada2.filter_and_trim(fwd, filt_f, rev, filt_r, **filt_kw, n_jobs=nt),
        )
    for nt in (1, 4):
        run(f"learn_errors@{nt}", lambda nt=nt: dada2.learn_errors(filt_f, n_jobs=nt))
        run(f"dada@{nt}", lambda nt=nt: dada2.DADA(err=ef, n_jobs=nt).fit(drp_f))
        run(f"dada_pool@{nt}", lambda nt=nt: dada2.DADA(err=ef, pool=True, n_jobs=nt).fit(drp_f))
        run(
            f"dada_pseudo@{nt}",
            lambda nt=nt: dada2.DADA(err=ef, pool="pseudo", n_jobs=nt).fit(drp_f),
        )
        for method in ("consensus", "pooled", "per-sample"):
            run(
                f"remove_bimera_{method.replace('-', '_')}@{nt}",
                lambda nt=nt, method=method: dada2.remove_bimera_denovo(
                    st_chim, method=method, n_jobs=nt
                ),
            )
        run(
            f"assign_taxonomy@{nt}",
            lambda nt=nt: dada2.assign_taxonomy(qfull, train, random_state=1, n_jobs=nt),
        )
        # single-threaded in R: the "@4" rows are pydada2-only
        run(
            f"merge_pairs@{nt}",
            lambda nt=nt: dada2.merge_pairs(dd_f, drp_f, dd_r, drp_r, n_jobs=nt),
        )
        run(
            f"collapse_no_mismatch@{nt}",
            lambda nt=nt: dada2.collapse_no_mismatch(st_chim, n_jobs=nt),
        )
        run(
            f"assign_species@{nt}",
            lambda nt=nt: dada2.assign_species(qspec, species_ref, n_jobs=nt),
        )
        run(
            f"remove_primers@{nt}",
            lambda nt=nt: dada2.remove_primers(
                example_path("samPBprimers.fastq.gz"), tmp / "p.fq.gz", f27, r1492rc, n_jobs=nt
            ),
        )
        run(
            f"remove_primers_indels@{nt}",
            lambda nt=nt: dada2.remove_primers(
                example_path("samPBprimers.fastq.gz"), tmp / "p2.fq.gz", f27, r1492rc,
                allow_indels=True, n_jobs=nt,
            ),
        )  # fmt: skip
        run(f"seq_complexity@{nt}", lambda nt=nt: dada2.seq_complexity(sam1_seqs, n_jobs=nt))
    run("derep_fastq@1", lambda: dada2.derep_fastq(filt_f))
    run(
        "make_sequence_table@1",
        lambda: dada2.make_sequence_table(dict(zip(["sam1", "sam2"], mergers, strict=True))),
    )
    res["meta"]["loops"] = LOOPS  # type: ignore[index]
    out.write_text(json.dumps(res, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
