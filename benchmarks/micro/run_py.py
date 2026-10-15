"""Micro-benchmarks of pydada2 on the bundled example data (spec §9.1).

Usage: uv run python benchmarks/micro/run_py.py [repeats] [outfile]

Stage names and inputs mirror benchmarks/micro/run_r.R. Each stage runs once as warm-up, then
``repeats`` times; the median wall time is reported.
"""

from __future__ import annotations

import json
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


def bench(fn: Callable[[], object], repeats: int) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fn()
        times = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            fn()
            times.append(time.perf_counter() - t0)
    return statistics.median(times)


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
    ten = [s for _, s in read_fasta(example_path("ten_16s.100.fa.gz"))]
    qfull = ten[::10]
    qspec = [s for s in qfull if set(s) <= set("ACGT")]
    f27, r1492rc = "AGRGTTYGATYMTGGCTCAG", dada2.reverse_complement("RGYTACCTTGTTACGACTT")
    train = example_path("example_train_set.fa.gz")
    sam1_seqs = dada2.get_sequences(example_path("sam1F.fastq.gz"))

    res: dict[str, object] = {
        "meta": {
            "python": platform.python_version(),
            "pydada2": dada2.__version__,
            "repeats": repeats,
            "cores": os.cpu_count(),
        }
    }
    for nt in (1, 4):

        def k(name: str, nt: int = nt) -> str:
            return f"{name}@{nt}"

        res[k("filter_and_trim")] = bench(
            lambda nt=nt: dada2.filter_and_trim(fwd, filt_f, rev, filt_r, **filt_kw, n_jobs=nt),
            repeats,
        )
        res[k("learn_errors")] = bench(lambda nt=nt: dada2.learn_errors(filt_f, n_jobs=nt), repeats)
        res[k("dada")] = bench(lambda nt=nt: dada2.DADA(err=ef, n_jobs=nt).fit(drp_f), repeats)
        res[k("dada_pool")] = bench(
            lambda nt=nt: dada2.DADA(err=ef, pool=True, n_jobs=nt).fit(drp_f), repeats
        )
        res[k("dada_pseudo")] = bench(
            lambda nt=nt: dada2.DADA(err=ef, pool="pseudo", n_jobs=nt).fit(drp_f), repeats
        )
        res[k("remove_bimera_consensus")] = bench(
            lambda nt=nt: dada2.remove_bimera_denovo(st, method="consensus", n_jobs=nt), repeats
        )
        res[k("remove_bimera_pooled")] = bench(
            lambda nt=nt: dada2.remove_bimera_denovo(st, method="pooled", n_jobs=nt), repeats
        )
        res[k("assign_taxonomy")] = bench(
            lambda nt=nt: dada2.assign_taxonomy(qfull, train, random_state=1, n_jobs=nt), repeats
        )
    res["derep_fastq@1"] = bench(lambda: dada2.derep_fastq(filt_f), repeats)
    res["merge_pairs@1"] = bench(lambda: dada2.merge_pairs(dd_f, drp_f, dd_r, drp_r), repeats)
    res["make_sequence_table@1"] = bench(
        lambda: dada2.make_sequence_table(dict(zip(["sam1", "sam2"], mergers, strict=True))),
        repeats,
    )
    res["collapse_no_mismatch@1"] = bench(lambda: dada2.collapse_no_mismatch(st), repeats)
    res["assign_species@1"] = bench(
        lambda: dada2.assign_species(qspec, example_path("example_species_assignment.fa.gz")),
        repeats,
    )
    res["remove_primers@1"] = bench(
        lambda: dada2.remove_primers(
            example_path("samPBprimers.fastq.gz"), tmp / "p.fq.gz", f27, r1492rc
        ),
        repeats,
    )
    res["remove_primers_indels@1"] = bench(
        lambda: dada2.remove_primers(
            example_path("samPBprimers.fastq.gz"),
            tmp / "p2.fq.gz",
            f27,
            r1492rc,
            allow_indels=True,
        ),
        repeats,
    )
    res["seq_complexity@1"] = bench(lambda: dada2.seq_complexity(sam1_seqs), repeats)
    out.write_text(json.dumps(res, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
