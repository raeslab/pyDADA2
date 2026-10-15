"""Run the R and Python MiSeq SOP pipelines repeatedly (spec §9.2, phase 11).

Usage: uv run python benchmarks/dataset/run.py [--tool R|py] [repeats] [threads ...]
       (defaults: both tools, 3 repeats at 1 and 8 threads)

Each run is a fresh process writing benchmarks/dataset/results/<tool>_t<threads>_r<rep>/.
Runs are interleaved (R, Python, R, Python, ...) so slow drifts of the machine affect both
tools alike. Peak RSS is the ``ru_maxrss`` reported by ``wait4`` for the run's process tree
(which includes forked workers, e.g. R's ``mclapply`` in ``filterAndTrim``).
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


def run(tool: str, threads: int, rep: int) -> None:
    out = RESULTS / f"{tool}_t{threads}_r{rep}"
    script = HERE / ("pipeline.R" if tool == "R" else "pipeline.py")
    cmd = ["Rscript", str(script)] if tool == "R" else [sys.executable, str(script)]
    cmd += [str(threads), str(out), str(HERE / "data")]
    t0 = time.perf_counter()
    with open(RESULTS / f"{out.name}.log", "w") as log:
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
    if args[:1] == ["--tool"]:
        tools, args = (args[1],), args[2:]
    repeats = int(args[0]) if args else 3
    threads = [int(t) for t in args[1:]] or [1, 8]
    RESULTS.mkdir(parents=True, exist_ok=True)
    for rep in range(1, repeats + 1):
        for t in threads:
            for tool in tools:
                run(tool, t, rep)


if __name__ == "__main__":
    main()
