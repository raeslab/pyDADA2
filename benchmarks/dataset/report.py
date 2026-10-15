"""Compare and summarise the MiSeq SOP benchmark runs (spec §9.2, phase 11).

Usage: uv run python benchmarks/dataset/report.py [results_dir] [report.md]

Reads every <tool>_t<threads>_r<rep>/results.json (+ run.json) written by run.py, checks the
outputs of pydada2 against R dada2 (exact for everything except taxonomy, which follows the
statistical criterion of spec §8.3), and writes a Markdown report with median timings.
"""

from __future__ import annotations

import json
import os
import platform
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
STAGES = [
    "filter_and_trim", "learn_errors_F", "learn_errors_R", "derep_F", "derep_R", "dada_F",
    "dada_R", "merge_pairs", "make_sequence_table", "remove_bimera_denovo", "assign_taxonomy",
    "add_species", "total",
]  # fmt: skip
EXACT = [
    ("read tracking table", "track"),
    ("unique sequences per sample (derep)", "derep"),
    ("denoised ASVs + abundances per sample (F and R)", "dada"),
    ("merged pairs (sequence, abundance, F/R ASV index)", "mergers"),
    ("sequence table (ASVs, column order, counts)", "seqtab"),
    ("chimera calls", "bimera"),
]


def as_list(x: Any) -> list[Any]:
    """R's auto_unbox turns length-1 vectors into scalars; undo that."""
    return x if isinstance(x, list) else [x]


def load(run_dir: Path) -> dict[str, Any]:
    d = json.loads((run_dir / "results.json").read_text())
    d["run"] = json.loads((run_dir / "run.json").read_text())
    return d


def per_sample(d: dict[str, Any], key: str, fields: list[str]) -> dict[str, list[tuple[Any, ...]]]:
    return {s: list(zip(*(as_list(v[f]) for f in fields), strict=True)) for s, v in d[key].items()}


def exact_checks(a: dict[str, Any], b: dict[str, Any]) -> dict[str, bool]:
    return {
        "track": a["samples"] == b["samples"]
        and all(as_list(a["track"][k]) == as_list(b["track"][k]) for k in a["track"]),
        "derep": a["derep_F"] == b["derep_F"] and a["derep_R"] == b["derep_R"],
        "dada": all(
            per_sample(a, k, ["seqs", "abund"]) == per_sample(b, k, ["seqs", "abund"])
            for k in ("dada_F", "dada_R")
        ),
        "mergers": per_sample(a, "mergers", ["sequence", "abundance", "forward", "reverse"])
        == per_sample(b, "mergers", ["sequence", "abundance", "forward", "reverse"]),
        "seqtab": a["seqtab"] == b["seqtab"],
        "bimera": as_list(a["bimera"]) == as_list(b["bimera"]),
    }


def err_rtol(a: dict[str, Any], b: dict[str, Any]) -> float:
    worst = 0.0
    for k in ("err_F", "err_R"):
        x, y = np.array(a[k], dtype=float), np.array(b[k], dtype=float)
        worst = max(worst, float(np.max(np.abs(x - y) / np.maximum(np.abs(y), 1e-300))))
    return worst


def bound(p: Any) -> Any:
    """3 sigma of the difference of two 100-draw bootstrap estimates (+2 for ties), as in
    tests/parity/test_taxonomy.py."""
    p = np.clip(np.asarray(p, dtype=float) / 100.0, 0.0, 1.0)
    return 3 * np.sqrt(2 * 100 * p * (1 - p)) + 2


def taxonomy(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    """Compare b against reference a (a: R)."""
    assert a["taxa"]["seqs"] == b["taxa"]["seqs"] and a["taxa"]["levels"] == b["taxa"]["levels"]
    levels = a["taxa"]["levels"]
    n = len(a["taxa"]["seqs"])
    ta = np.array(a["taxa"]["values"], dtype=object).reshape(n, len(levels))
    tb = np.array(b["taxa"]["values"], dtype=object).reshape(n, len(levels))
    nb = len(a["boot"]["levels"])
    ba = np.array(a["boot"]["values"], dtype=float).reshape(n, nb)
    bb = np.array(b["boot"]["values"], dtype=float).reshape(n, nb)
    same = ta == tb
    # Disagreements at min_boot=50 are allowed only where a's bootstrap is within noise of 50
    # (first differing level of each sequence; deeper levels follow from it).
    unexplained, borderline = 0, []
    for i in range(n):
        for j in range(nb):
            if not same[i, j]:
                borderline.append(ba[i, j])
                if abs(ba[i, j] - 50) > bound(50.0):
                    unexplained += 1
                break
    diff = np.abs(ba - bb)
    return {
        "n": n,
        "levels": levels,
        "identical_per_level": same.sum(axis=0).tolist(),
        "identical_rows": int(same.all(axis=1).sum()),
        "boot_median": float(np.median(diff)),
        "boot_max": float(diff.max()),
        "boot_outside_3sigma": int((diff > bound((ba + bb) / 2)).sum()),
        "borderline": borderline,
        "unexplained": unexplained,
    }


def fmt_t(s: float) -> str:
    return f"{s / 60:.1f} min" if s >= 120 else f"{s:.1f} s" if s >= 1 else f"{s * 1e3:.0f} ms"


def cpu_name() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine()


def main() -> None:
    results = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "results"
    report = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "REPORT.md"
    runs: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for p in sorted(results.iterdir()):
        m = re.fullmatch(r"(R|py)_t(\d+)_r(\d+)", p.name)
        if m and (p / "run.json").exists():
            runs[(m[1], int(m[2]))].append(load(p))
    threads = sorted({t for _, t in runs})
    ref, py = runs[("R", threads[0])][0], runs[("py", threads[0])][0]

    lines: list[str] = []
    lines += [
        "# MiSeq SOP benchmark: pydada2 vs R dada2",
        "",
        "Generated by `benchmarks/dataset/report.py` from runs of `run.py` (spec §9.2, phase 11).",
        "Both tools run the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html)",
        "pipeline (`pipeline.R` / `pipeline.py`) with identical parameters, each run in a fresh",
        "process.",
        "",
        "* Data: MiSeq SOP (Kozich et al. 2013), 20 samples, 16S V4, 2x250 "
        f"({sum(as_list(ref['track']['input'])):,} read pairs); SILVA 138.2 DADA2 training "
        "set (to genus) and species set",
        f"* R: {ref['meta']['version']}",
        f"* Python: {py['meta']['version']}",
        f"* Machine: {cpu_name()}, {os.cpu_count()} logical CPUs, {platform.system()}",
        f"* Runs: {', '.join(f'{len(v)} x {k[0]} @ {k[1]} thr' for k, v in sorted(runs.items()))}",
        "",
        "## Results: identical",
        "",
        f"pydada2 vs R at {threads[0]} thread(s), first run of each:",
        "",
        "| output | identical |",
        "|---|---|",
    ]
    checks = exact_checks(ref, py)
    lines += [f"| {label} | {'yes' if checks[k] else '**NO**'} |" for label, k in EXACT]
    lines.append(
        f"| error rates (learn_errors F and R) | max rel. difference {err_rtol(py, ref):.1e} |"
    )

    # Determinism across repeats and thread counts, per tool.
    lines += ["", "Every run of each tool, at every thread count, gives the same results as that "
          "tool's first run:", ""]  # fmt: skip
    for tool in ("R", "py"):
        first = runs[(tool, threads[0])][0]
        others = [r for t in threads for r in runs[(tool, t)]]
        same = all(all(exact_checks(first, r).values()) for r in others)
        lines.append(f"* {'R' if tool == 'R' else 'pydada2'}: {'yes' if same else '**no**'} "
                 f"({len(others)} runs, exact outputs)")  # fmt: skip

    # Taxonomy: R vs pydada2 next to R vs R (upstream is not reproducible, see DEVIATIONS.md).
    tx = taxonomy(ref, py)
    rr = runs[("R", threads[0])]
    tx_rr = taxonomy(rr[0], rr[1]) if len(rr) > 1 else None
    lines += [
        "",
        "### Taxonomy",
        "",
        "Bootstrap confidences use different random numbers (DEVIATIONS.md #2), and R itself",
        "breaks ties with `std::random_device`, so R runs are not exactly reproducible either",
        "(on tied genera). The criterion is therefore statistical (spec §8.3).",
        "",
        f"| {tx['n']} non-chimeric ASVs | pydada2 vs R |" + (" R vs R |" if tx_rr else ""),
        "|---|---|" + ("---|" if tx_rr else ""),
    ]
    for j, lev in enumerate(tx["levels"]):
        row = f"| identical {lev} | {tx['identical_per_level'][j]} |"
        lines.append(row + (f" {tx_rr['identical_per_level'][j]} |" if tx_rr else ""))
    for label, k in [
        ("identical at all levels", "identical_rows"),
        ("median \\|Δ bootstrap\\|", "boot_median"),
        ("max \\|Δ bootstrap\\|", "boot_max"),
        ("entries outside 3 sigma binomial noise", "boot_outside_3sigma"),
        ("disagreements with R bootstrap far from 50", "unexplained"),
    ]:
        v = tx[k]
        row = f"| {label} | {v:g} |"
        lines.append(row + (f" {tx_rr[k]:g} |" if tx_rr else ""))
    if tx["borderline"]:
        b = tx["borderline"]
        lines += [
            "",
            f"The R bootstraps at the levels where pydada2 differs are {min(b):.0f} to "
            f"{max(b):.0f}, close to the `min_boot=50` cutoff: borderline calls.",
        ]

    # Timings.
    lines += [
        "",
        "## Speed",
        "",
        "Median wall time per stage over the repeats (ratio = R / pydada2, > 1 means pydada2 "
        "is faster).",
        "",
    ]
    head = "| stage |" + "".join(f" R {t} thr | pydada2 {t} thr | speed-up |" for t in threads)
    lines += [head, "|---|" + "---|---|---|" * len(threads)]

    def med(tool: str, t: int, f: Any) -> float:
        return statistics.median(f(r) for r in runs[(tool, t)])

    for s in STAGES:
        row = f"| {'**total**' if s == 'total' else s} |"
        for t in threads:
            r = med("R", t, lambda d, s=s: d["times"][s])
            p = med("py", t, lambda d, s=s: d["times"][s])
            row += f" {fmt_t(r)} | {fmt_t(p)} | {r / p:.1f}x |"
        lines.append(row)
    for label, f in [
        ("wall time incl. start-up", lambda d: d["run"]["wall_incl_startup"]),
    ]:
        row = f"| {label} |"
        for t in threads:
            r, p = med("R", t, f), med("py", t, f)
            row += f" {fmt_t(r)} | {fmt_t(p)} | {r / p:.1f}x |"
        lines.append(row)
    row = "| peak RSS |"
    for t in threads:
        r = med("R", t, lambda d: d["run"]["max_rss_kib"]) / 1024**2
        p = med("py", t, lambda d: d["run"]["max_rss_kib"]) / 1024**2
        row += f" {r:.2f} GiB | {p:.2f} GiB | {r / p:.1f}x |"
    lines.append(row)
    report.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    ok = all(checks.values()) and tx["unexplained"] == 0
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
