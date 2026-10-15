"""Run the R and Python MiSeq SOP pipelines repeatedly (spec §9.2, phase 11).

Usage: uv run python benchmarks/dataset/run.py [--tool R|py] [--fastq-dir NAME] [--results DIR]
                                                 [--python EXE]
                                                 [repeats] [threads ...]
       (defaults: both tools, MiSeq_SOP, benchmarks/dataset/results, 3 repeats at 1 and 8 threads)

Each run is a fresh process writing <results>/<tool>_t<threads>_r<rep>/. ``--fastq-dir`` names
the FASTQ sub-directory of benchmarks/dataset/data (e.g. Stability, see download.sh);
``--python`` the interpreter for the pydada2 runs (e.g. another venv with an older pydada2).
Runs are interleaved (R, Python, R, Python, ...) so slow drifts of the machine affect both
tools alike. Peak RSS is the ``ru_maxrss`` reported by ``wait4``: the peak resident set of the
*largest single process* in the run's process tree (the pipeline itself or one of its waited-for
children, e.g. a forked ``mclapply`` worker of R's ``filterAndTrim``), not the sum over
concurrently running processes. Threads share their process's memory, so for pydada2 (threads
only) it is the whole footprint; for R runs with forked workers the total can be higher.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def run(tool: str, threads: int, rep: int, fastq_dir: str, results: Path, python: str) -> None:
    out = results / f"{tool}_t{threads}_r{rep}"
    script = HERE / ("pipeline.R" if tool == "R" else "pipeline.py")
    cmd = ["Rscript", str(script)] if tool == "R" else [python, str(script)]
    cmd += [str(threads), str(out), str(HERE / "data"), fastq_dir]
    t0 = time.perf_counter()
    with open(results / f"{out.name}.log", "w") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
    wall = time.perf_counter() - t0
    if os.waitstatus_to_exitcode(status) != 0:
        raise SystemExit(f"{out.name} failed, see {out.name}.log")
    meta = {"wall_incl_startup": wall, "max_rss_kib": usage.ru_maxrss}
    (out / "run.json").write_text(json.dumps(meta))
    print(f"{out.name}: {wall:.1f} s, peak RSS {usage.ru_maxrss / 1024**2:.2f} GiB", flush=True)


def main() -> None:
    args = sys.argv[1:]
    tools: tuple[str, ...] = ("R", "py")
    fastq_dir, results, python = "MiSeq_SOP", RESULTS, sys.executable
    while args[:1] and args[0].startswith("--"):
        opt, val, args = args[0], args[1], args[2:]
        if opt == "--tool":
            tools = (val,)
        elif opt == "--fastq-dir":
            fastq_dir = val
        elif opt == "--results":
            results = Path(val).resolve()
        elif opt == "--python":
            python = val
        else:
            raise SystemExit(f"unknown option {opt}")
    repeats = int(args[0]) if args else 3
    threads = [int(t) for t in args[1:]] or [1, 8]
    results.mkdir(parents=True, exist_ok=True)
    for rep in range(1, repeats + 1):
        for t in threads:
            for tool in tools:
                run(tool, t, rep, fastq_dir, results, python)


if __name__ == "__main__":
    main()
