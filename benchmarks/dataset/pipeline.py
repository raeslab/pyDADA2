"""DADA2 tutorial pipeline on the MiSeq SOP data with pydada2 (spec §9.2, phase 11).

Usage: python benchmarks/dataset/pipeline.py <threads> <outdir> [datadir] [fastq_subdir]

fastq_subdir defaults to MiSeq_SOP; files ending in _R1_001.fastq(.gz) / _R2_001.fastq(.gz) are
used. Samples with no reads left after filtering are dropped (neither tool writes their files).

Mirrors benchmarks/dataset/pipeline.R stage by stage, with the same parameters. Writes
<outdir>/results.json in the same layout (indices 1-based, matrices row-major, NA as null).
"""

from __future__ import annotations

import json
import platform
import re
import sys
import time
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

import pydada2 as dada2


def main() -> None:
    threads = int(sys.argv[1])
    out = Path(sys.argv[2])
    data = Path(sys.argv[3] if len(sys.argv) > 3 else "benchmarks/dataset/data")
    out.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")

    times: dict[str, float] = {}

    def stage(name: str, fn: Callable[[], Any]) -> Any:
        t0 = time.perf_counter()
        val = fn()
        times[name] = time.perf_counter() - t0
        return val

    sop = data / (sys.argv[4] if len(sys.argv) > 4 else "MiSeq_SOP")
    fn_fs = sorted(p for p in sop.iterdir() if re.search(r"_R1_001\.fastq(\.gz)?$", p.name))
    fn_rs = sorted(p for p in sop.iterdir() if re.search(r"_R2_001\.fastq(\.gz)?$", p.name))
    samples = [p.name.split("_")[0] for p in fn_fs]
    filt_fs = [out / "filtered" / f"{s}_F_filt.fastq.gz" for s in samples]
    filt_rs = [out / "filtered" / f"{s}_R_filt.fastq.gz" for s in samples]

    t_total = time.perf_counter()
    stats = stage(
        "filter_and_trim",
        lambda: dada2.filter_and_trim(
            fn_fs, filt_fs, fn_rs, filt_rs, trunc_len=(240, 160), max_n=0, max_ee=(2, 2),
            trunc_q=2, rm_phix=True, compress=True, n_jobs=threads,
        ),
    )  # fmt: skip
    # Neither tool writes the files of samples without reads left: drop those samples.
    keep = (stats["reads_out"] > 0).tolist()
    samples = [s for s, k in zip(samples, keep, strict=True) if k]
    filt_fs = [f for f, k in zip(filt_fs, keep, strict=True) if k]
    filt_rs = [f for f, k in zip(filt_rs, keep, strict=True) if k]
    err_f = stage("learn_errors_F", lambda: dada2.learn_errors(filt_fs, n_jobs=threads))
    err_r = stage("learn_errors_R", lambda: dada2.learn_errors(filt_rs, n_jobs=threads))
    drp_f = stage("derep_F", lambda: dict(zip(samples, dada2.derep_fastq(filt_fs), strict=True)))
    drp_r = stage("derep_R", lambda: dict(zip(samples, dada2.derep_fastq(filt_rs), strict=True)))
    dd_f = stage("dada_F", lambda: dada2.DADA(err=err_f.err_, n_jobs=threads).fit(drp_f))
    dd_r = stage("dada_R", lambda: dada2.DADA(err=err_r.err_, n_jobs=threads).fit(drp_r))
    mergers = stage(
        "merge_pairs",
        lambda: dada2.merge_pairs(
            dd_f, list(drp_f.values()), dd_r, list(drp_r.values()), n_jobs=threads
        ),
    )
    seqtab = stage(
        "make_sequence_table",
        lambda: dada2.make_sequence_table(dict(zip(samples, mergers, strict=True))),
    )
    remover = stage(
        "remove_bimera_denovo",
        lambda: dada2.BimeraRemover(method="consensus", n_jobs=threads).fit(seqtab),
    )
    nochim = remover.transform(seqtab)
    tax, boot = stage(
        "assign_taxonomy",
        lambda: dada2.assign_taxonomy(
            list(nochim.columns), data / "silva_nr99_v138.2_toGenus_trainset.fa.gz",
            output_bootstraps=True, random_state=100, n_jobs=threads,
        ),
    )  # fmt: skip
    species = stage(
        "add_species",
        lambda: dada2.add_species(tax, data / "silva_v138.2_assignSpecies.fa.gz", n_jobs=threads),
    )
    times["total"] = time.perf_counter() - t_total

    def uniq(dd: dada2.DADA) -> dict[str, Any]:
        return {
            s: {"seqs": r.sequence, "abund": r.denoised.tolist()}
            for s, r in zip(samples, dd.results_, strict=True)
        }

    def values(df: Any) -> list[Any]:
        return [None if v is None or v != v else v for v in df.to_numpy(object).ravel().tolist()]

    res = {
        "meta": {
            "tool": "pydada2",
            "version": f"Python {platform.python_version()} / pydada2 {dada2.__version__}",
            "threads": threads,
        },
        "times": times,
        "samples": samples,
        "track": {
            "input": stats["reads_in"].tolist(),
            "filtered": stats["reads_out"].tolist(),
            "denoisedF": [int(r.denoised.sum()) for r in dd_f.results_],
            "denoisedR": [int(r.denoised.sum()) for r in dd_r.results_],
            "merged": [int(m["abundance"][m["accept"]].sum()) for m in mergers],
            "nonchim": nochim.sum(axis=1).astype(int).tolist(),
        },
        "err_F": np.asarray(dada2.get_errors(err_f)).ravel().tolist(),
        "err_R": np.asarray(dada2.get_errors(err_r)).ravel().tolist(),
        "derep_F": [len(d) for d in drp_f.values()],
        "derep_R": [len(d) for d in drp_r.values()],
        "dada_F": uniq(dd_f),
        "dada_R": uniq(dd_r),
        "mergers": {
            s: {
                "sequence": m["sequence"].tolist(),
                "abundance": m["abundance"].tolist(),
                "forward": (m["forward"] + 1).tolist(),  # R is 1-based
                "reverse": (m["reverse"] + 1).tolist(),
            }
            for s, m in zip(samples, mergers, strict=True)
        },
        "seqtab": {
            "seqs": [str(c) for c in seqtab.columns],
            "counts": seqtab.to_numpy().ravel().tolist(),
        },
        "bimera": remover.is_bimera_.tolist(),
        "taxa": {
            "seqs": [str(i) for i in species.index],
            "levels": list(species.columns),
            "values": values(species),
        },
        "boot": {"levels": list(boot.columns), "values": boot.to_numpy().ravel().tolist()},
    }
    (out / "results.json").write_text(json.dumps(res))


if __name__ == "__main__":
    main()
