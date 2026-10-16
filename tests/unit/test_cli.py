"""Command-line interface (pydada2.cli) and its file formats (pydada2.formats)."""

import gzip
import pickle
import subprocess
import sys

import click
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import pydada2 as dada2
from pydada2.cli import EXIT_NO_DATA, GreedyOption, cli, main, sample_name_from_path
from pydada2.datasets import err_balanced_f, example_path
from pydada2.formats import _empty_dada

NAMES = ("sam1", "sam2")
TRAIN = example_path("example_train_set.fa.gz")
SPECIES = example_path("example_species_assignment.fa.gz")


def run(*args, code=0):
    res = CliRunner().invoke(cli, [str(a) for a in args])
    assert res.exit_code == code, (res.output, res.stderr, res.exception)
    return res


def tsv(text):
    lines = text.strip("\n").split("\n")
    return [dict(zip(lines[0].split("\t"), ln.split("\t"), strict=True)) for ln in lines[1:]]


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    """The whole CLI workflow on the example data (paths and stdout of every step)."""
    d = tmp_path_factory.mktemp("cli")
    out = {"dir": d, "stdout": {}}
    for n in NAMES:
        res = run(
            "trim-and-filter",
            "--in-fwd", example_path(f"{n}F.fastq.gz"), "--in-rev", example_path(f"{n}R.fastq.gz"),
            "--out-fwd", d / "filt" / f"{n}_F.fq.gz", "--out-rev", d / "filt" / f"{n}_R.fq.gz",
            "--trunc-len", 240, 160, "--max-ee", 2, "--sample-name", n,
            "--stats", d / "stats" / f"{n}.tsv.gz",
        )  # fmt: skip
        out["stdout"][f"filter_{n}"] = res.stdout
    filt_f = [d / "filt" / f"{n}_F.fq.gz" for n in NAMES]
    filt_r = [d / "filt" / f"{n}_R.fq.gz" for n in NAMES]
    res = run(
        "learn-errors", "--in-fwd", *filt_f, "--in-rev", *filt_r, "--out", d / "errors.pkl",
        "--out-tsv", d / "errors.tsv.gz", "--plot", d / "errors.png", "--threads", 2,
    )  # fmt: skip
    out["stdout"]["learn"] = res.stdout
    for n in NAMES:
        res = run(
            "infer-asv", "--in-fwd", d / "filt" / f"{n}_F.fq.gz",
            "--in-rev", d / "filt" / f"{n}_R.fq.gz", "--errors", d / "errors.pkl",
            "--out-fwd", d / "asv" / f"{n}_F.dada.json.gz",
            "--out-rev", d / "asv" / f"{n}_R.dada.json.gz",
        )  # fmt: skip
        out["stdout"][f"infer_{n}"] = res.stdout
    # reverse files listed in the opposite order: pairing is by the stored sample name
    res = run(
        "merge",
        "--in-fwd", *[d / "asv" / f"{n}_F.dada.json.gz" for n in NAMES],
        "--in-rev", *[d / "asv" / f"{n}_R.dada.json.gz" for n in reversed(NAMES)],
        "--out-seqtab", d / "seqtab.tsv", "--out-seqtab-nochim", d / "nochim.tsv",
        "--out-fasta", d / "asvs.fa", "--out-mergers", d / "mergers.tsv", "-v",
    )  # fmt: skip
    out["stdout"]["merge"] = res.stdout
    res = run(
        "assign-taxonomy", "--in", d / "nochim.tsv", "--ref", TRAIN,
        "--out-taxa", d / "taxa.tsv", "--out-bootstrap", d / "boot.tsv",
        "--species-ref", SPECIES, "--save-model", d / "model.pkl",
    )  # fmt: skip
    out["stdout"]["taxonomy"] = res.stdout
    return out


@pytest.fixture(scope="module")
def api(pipeline):
    """The same workflow through the Python API."""
    d = pipeline["dir"]
    fwd = [example_path(f"{n}F.fastq.gz") for n in NAMES]
    rev = [example_path(f"{n}R.fastq.gz") for n in NAMES]
    filt_f = [d / "api" / f"{n}_F.fq.gz" for n in NAMES]
    filt_r = [d / "api" / f"{n}_R.fq.gz" for n in NAMES]
    stats = dada2.filter_and_trim(fwd, filt_f, rev, filt_r, trunc_len=(240, 160), max_ee=2)
    ef, er = dada2.learn_errors(filt_f), dada2.learn_errors(filt_r)
    drp_f, drp_r = dada2.derep_fastq(filt_f), dada2.derep_fastq(filt_r)
    dd_f = dada2.DADA(err=ef.err_).fit(drp_f)
    dd_r = dada2.DADA(err=er.err_).fit(drp_r)
    mergers = dada2.merge_pairs(dd_f, drp_f, dd_r, drp_r)
    seqtab = dada2.make_sequence_table(dict(zip(NAMES, mergers, strict=True)))
    nochim = dada2.BimeraRemover().fit_transform(seqtab)
    clf = dada2.TaxonomyClassifier(random_state=0).fit(TRAIN)
    taxa, boot = clf.predict(list(nochim.columns), return_bootstrap=True)
    taxa = dada2.add_species(taxa, SPECIES)
    return dict(stats=stats, ef=ef, er=er, dd_f=dd_f, dd_r=dd_r, drp_f=drp_f, drp_r=drp_r,
                mergers=mergers, seqtab=seqtab, nochim=nochim, taxa=taxa, boot=boot)  # fmt: skip


# ------------------------------------------------------------------ parity with the API ---
def test_trim_and_filter_matches_api(pipeline, api):
    d = pipeline["dir"]
    for i, n in enumerate(NAMES):
        row = tsv(pipeline["stdout"][f"filter_{n}"])
        assert row == [
            {"sample": n, "reads_in": str(api["stats"]["reads_in"].iloc[i]),
             "reads_out": str(api["stats"]["reads_out"].iloc[i])}
        ]  # fmt: skip
        with gzip.open(d / "stats" / f"{n}.tsv.gz", "rt") as fh:
            assert fh.read() == pipeline["stdout"][f"filter_{n}"]
        for m in "FR":
            cli_reads = next(dada2.io.read_fastq_chunks(d / "filt" / f"{n}_{m}.fq.gz"))
            api_reads = next(dada2.io.read_fastq_chunks(d / "api" / f"{n}_{m}.fq.gz"))
            assert (cli_reads.ids, cli_reads.seqs, cli_reads.quals) == (
                api_reads.ids, api_reads.seqs, api_reads.quals,
            )  # fmt: skip


def test_learn_errors_matches_api(pipeline, api):
    d = pipeline["dir"]
    with open(d / "errors.pkl", "rb") as fh:
        err = pickle.load(fh)
    np.testing.assert_array_equal(err["err_f"], api["ef"].err_)
    np.testing.assert_array_equal(err["err_r"], api["er"].err_)
    rows = tsv(pipeline["stdout"]["learn"])
    assert [r["mate"] for r in rows] == ["fwd", "rev"]
    assert rows[0]["n_reads"] == str(api["ef"].n_reads_)
    assert rows[1]["n_iter"] == str(api["er"].n_iter_)
    assert rows[0]["converged"] == ("true" if api["ef"].converged_ else "false")
    long = pd.read_csv(d / "errors.tsv.gz", sep="\t", float_precision="round_trip")
    assert list(long.columns) == ["mate", "transition", "qual", "rate"]
    sub = long[(long["mate"] == "rev") & (long["transition"] == "A2C")]
    np.testing.assert_array_equal(sub["rate"].to_numpy(), api["er"].err_[1])
    assert (d / "errors.png").stat().st_size > 0


def test_infer_asv_matches_api(pipeline, api):
    d = pipeline["dir"]
    for i, n in enumerate(NAMES):
        for mate, dd, drp in (("F", api["dd_f"], api["drp_f"]), ("R", api["dd_r"], api["drp_r"])):
            s = dada2.read_dada_json(d / "asv" / f"{n}_{mate}.dada.json.gz")
            assert (s.sample, s.mate) == (n, "fwd" if mate == "F" else "rev")
            assert s.source == f"{n}_{mate}.fq.gz"
            assert s.pydada2_version == dada2.__version__
            pd.testing.assert_frame_equal(s.result.clusters, dd.results_[i].clusters)
            pd.testing.assert_frame_equal(s.result.birth_subs, dd.results_[i].birth_subs)
            np.testing.assert_array_equal(s.result.map, dd.results_[i].map)
            np.testing.assert_array_equal(s.result.pval, dd.results_[i].pval)
            np.testing.assert_array_equal(s.result.trans, dd.results_[i].trans)
            assert s.result.params == dd.results_[i].params
            pd.testing.assert_series_equal(s.derep.uniques, drp[i].uniques)
            np.testing.assert_array_equal(s.derep.map, drp[i].map)
        row = tsv(pipeline["stdout"][f"infer_{n}"])[0]
        assert row["denoised_fwd"] == str(int(api["dd_f"].results_[i].denoised.sum()))
        assert row["uniques_rev"] == str(len(api["drp_r"][i]))
        assert row["reads"] == str(len(api["drp_f"][i].map))


def test_merge_matches_api(pipeline, api):
    d = pipeline["dir"]
    seqtab, ids = dada2.read_seqtab(d / "seqtab.tsv", return_ids=True)
    nochim, ids_nc = dada2.read_seqtab(d / "nochim.tsv", return_ids=True)
    pd.testing.assert_frame_equal(seqtab, api["seqtab"].rename(index=dict(enumerate(NAMES))),
                                  check_index_type=False, check_column_type=False)  # fmt: skip
    assert list(nochim.columns) == list(api["nochim"].columns)
    assert list(ids) == [f"ASV{i + 1}" for i in range(seqtab.shape[1])]
    assert all(ids[s] == ids_nc[s] for s in nochim.columns)  # ids are kept
    fasta = dada2.io.read_fasta(d / "asvs.fa")
    assert fasta == [(ids_nc[s], s) for s in nochim.columns]
    mergers = pd.read_csv(d / "mergers.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(mergers.columns[:3]) == ["sample", "asv_id", "sequence"]
    assert set(mergers["accept"]) == {"true"}
    rows = tsv(pipeline["stdout"]["merge"])
    assert [r["sample"] for r in rows] == list(NAMES)
    for r, m in zip(rows, api["mergers"], strict=True):
        assert int(r["merged"]) == int(m["abundance"].sum())


def test_assign_taxonomy_matches_api(pipeline, api):
    d = pipeline["dir"]
    taxa = pd.read_csv(d / "taxa.tsv", sep="\t", dtype=str, keep_default_na=False)
    nochim = pd.read_csv(d / "nochim.tsv", sep="\t", dtype=str)
    assert taxa["asv_id"].tolist() == nochim["asv_id"].tolist()
    assert taxa["sequence"].tolist() == nochim["sequence"].tolist()
    expected = api["taxa"].fillna("").to_numpy()
    np.testing.assert_array_equal(taxa.iloc[:, 2:].to_numpy(), expected)
    boot = pd.read_csv(d / "boot.tsv", sep="\t")
    np.testing.assert_array_equal(boot.iloc[:, 2:].to_numpy(), api["boot"].to_numpy())
    rows = tsv(pipeline["stdout"]["taxonomy"])
    assert [r["rank"] for r in rows] == list(api["taxa"].columns)
    assert rows[0]["assigned"] == str(int(api["taxa"].iloc[:, 0].notna().sum()))
    # a saved model gives the same result (settings come from the command line)
    run("assign-taxonomy", "--in", d / "nochim.tsv", "--ref", d / "model.pkl",
        "--out-taxa", d / "taxa2.tsv", "--species-ref", SPECIES)  # fmt: skip
    assert (d / "taxa2.tsv").read_text() == (d / "taxa.tsv").read_text()


def test_python_m_and_console_script(pipeline):
    res = subprocess.run([sys.executable, "-m", "pydada2", "--help"], capture_output=True,
                         text=True, check=True)  # fmt: skip
    assert "trim-and-filter" in res.stdout
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0


# ------------------------------------------------------------------------- parsing ---
def test_greedy_option_values():
    @click.command()
    @click.option("--files", cls=GreedyOption)
    @click.option(
        "--pair",
        cls=GreedyOption,
        type=float,
        max_values=2,
        accept=lambda t: t.lstrip("-").replace(".", "").isdigit() or "inf" in t,
    )
    @click.option("--flag/--no-flag")
    def cmd(files, pair, flag):
        click.echo(f"{files}|{pair}|{flag}")  # fmt: skip

    def out(*args):
        res = CliRunner().invoke(cmd, list(args))
        assert res.exit_code == 0, res.output
        return res.output.strip()

    assert out("--files", "a", "b", "c", "--flag") == "('a', 'b', 'c')|()|True"
    assert out("--files=a", "b", "--files", "c") == "('a', 'b', 'c')|()|False"
    assert out("--pair", "-inf", "-4", "--files", "x") == "('x',)|(-inf, -4.0)|False"
    assert out("--pair", "1", "2") == "()|(1.0, 2.0)|False"
    assert CliRunner().invoke(cmd, ["--pair", "1", "2", "3"]).exit_code == 2  # 3 is not an option


@pytest.mark.parametrize(
    ("path", "name"),
    [
        ("x/sam1_R1.fastq.gz", "sam1"),
        ("sam1_1.fq", "sam1"),
        ("sam1_F.fq.gz", "sam1"),
        ("sam1_F.dada.json.gz", "sam1"),
        ("sam1F.fastq.gz", "sam1F"),
        ("_F.fq", "_F"),
        ("sample.txt", "sample.txt"),
    ],
)
def test_sample_name_from_path(path, name):
    assert sample_name_from_path(path) == name


def test_per_mate_values_and_inf(tmp_path):
    out = tmp_path / "o"
    base = ["trim-and-filter", "--in-fwd", example_path("sam1F.fastq.gz"),
            "--in-rev", example_path("sam1R.fastq.gz"), "--out-fwd", out / "F.fq",
            "--out-rev", out / "R.fq"]  # fmt: skip
    res = run(*base, "--min-len", "-inf", "--max-n", "inf", "--min-q", "-inf", 1,
              "--max-len", "inf", 300, "--trunc-len", 240, 160, "--rm-lowcomplex", 0,
              "--keep-phix")  # fmt: skip
    expected = dada2.filter_and_trim(
        example_path("sam1F.fastq.gz"), tmp_path / "F.fq", example_path("sam1R.fastq.gz"),
        tmp_path / "R.fq", compress=False, min_len=-np.inf, max_n=np.inf,
        min_q=(-np.inf, 1), max_len=(np.inf, 300), trunc_len=(240, 160), rm_phix=False,
    )  # fmt: skip
    assert tsv(res.stdout)[0]["reads_out"] == str(expected["reads_out"].iloc[0])
    assert (out / "F.fq").read_bytes() == (tmp_path / "F.fq").read_bytes()  # not gzipped
    assert tsv(res.stdout)[0]["sample"] == "sam1F"
    res = run(*base, "--trunc-len", 240, 160, 100, code=2)
    assert "extra argument (100)" in res.stderr
    res = run(*base, "--trunc-len", 240, 160, "--trunc-len", 100, code=2)
    assert "two values" in res.stderr
    res = run(*base, "--min-len", "abc", code=2)
    assert "inf" in res.stderr
    res = run(*base[:-2], "--out-rev", out / "R.fq.gz", code=2)
    assert ".gz" in res.stderr
    res = run(*base[:-2], "--out-rev", example_path("sam1R.fastq.gz"), code=2)
    assert "also an input" in res.stderr
    res = run(*base[:-2], "--out-rev", out / "F.fq", code=2)
    assert "distinct" in res.stderr
    run(*base[:-2], code=2)  # missing required option


# ---------------------------------------------------------------------- exit codes ---
def test_trim_and_filter_no_reads(tmp_path):
    res = run(
        "trim-and-filter", "--in-fwd", example_path("sam1F.fastq.gz"),
        "--in-rev", example_path("sam1R.fastq.gz"), "--out-fwd", tmp_path / "F.fq.gz",
        "--out-rev", tmp_path / "R.fq.gz", "--trunc-len", 1000, "--stats", tmp_path / "s.tsv",
        code=EXIT_NO_DATA,
    )  # fmt: skip
    assert tsv(res.stdout) == [{"sample": "sam1F", "reads_in": "1500", "reads_out": "0"}]
    assert (tmp_path / "s.tsv").read_text() == res.stdout
    assert not (tmp_path / "F.fq.gz").exists()
    assert "no reads passed" in res.stderr


def test_runtime_error_exit_1(tmp_path):
    args = ["trim-and-filter", "--in-fwd", example_path("sam1F.fastq.gz"),
            "--in-rev", example_path("sam1R.fastq.gz"), "--out-fwd", tmp_path / "F.fq.gz",
            "--out-rev", tmp_path / "R.fq.gz", "--orient-fwd", "NNNN"]  # fmt: skip
    res = run(*args, code=1)
    assert res.stderr.startswith("Error: Non-ACGT")
    res = run(*args, "-vv", code=1)
    assert "Traceback" in res.stderr


def test_infer_asv_empty_input(tmp_path, pipeline):
    for m in "FR":
        (tmp_path / f"e_{m}.fq").write_text("")
    res = run(
        "infer-asv", "--in-fwd", tmp_path / "e_F.fq", "--in-rev", tmp_path / "e_R.fq",
        "--errors", pipeline["dir"] / "errors.pkl", "--out-fwd", tmp_path / "e_F.dada.json",
        "--out-rev", tmp_path / "e_R.dada.json", "--band-size", 32, code=EXIT_NO_DATA,
    )  # fmt: skip
    assert tsv(res.stdout)[0] == {"sample": "e", "reads": "0", "uniques_fwd": "0",
                                  "uniques_rev": "0", "denoised_fwd": "0",
                                  "denoised_rev": "0"}  # fmt: skip
    s = dada2.read_dada_json(tmp_path / "e_F.dada.json")
    assert len(s.result.clusters) == 0 and len(s.derep) == 0
    assert s.result.params.band_size == 32
    # merging only empty samples: tables are written, exit code 3
    res = run(
        "merge", "--in-fwd", tmp_path / "e_F.dada.json", "--in-rev", tmp_path / "e_R.dada.json",
        "--out-seqtab", tmp_path / "st.tsv", "--out-seqtab-nochim", tmp_path / "nc.tsv",
        "--out-mergers", tmp_path / "m.tsv", code=EXIT_NO_DATA,
    )  # fmt: skip
    assert tsv(res.stdout) == [{"sample": "e", "merged": "0", "nonchim": "0"}]
    assert (tmp_path / "nc.tsv").read_text() == "asv_id\tsequence\te\n"
    # a table without ASVs: header-only taxonomy, exit code 3
    res = run("assign-taxonomy", "--in", tmp_path / "nc.tsv", "--ref", TRAIN,
              "--out-taxa", tmp_path / "t.tsv", "--out-bootstrap", tmp_path / "b.tsv",
              "--species-ref", SPECIES, code=EXIT_NO_DATA)  # fmt: skip
    assert (tmp_path / "t.tsv").read_text().split("\n")[0].endswith("Genus\tSpecies")
    assert all(r["assigned"] == "0" for r in tsv(res.stdout))


def test_infer_asv_options(tmp_path, pipeline):
    d = pipeline["dir"]
    with gzip.open(tmp_path / "err.pkl.gz", "wb") as fh:  # a gzipped pickle is accepted too
        pickle.dump({"err_f": err_balanced_f(), "err_r": err_balanced_f()}, fh)
    seqs = [r[1] for r in dada2.io.read_fasta(d / "asvs.fa")][:1]
    dada2.io.write_fasta([("p", seqs[0][:240])], tmp_path / "priors.fa")
    res = run(
        "infer-asv", "--in-fwd", d / "filt" / "sam1_F.fq.gz",
        "--in-rev", d / "filt" / "sam2_R.fq.gz",
        "--errors", tmp_path / "err.pkl.gz", "--out-fwd", tmp_path / "F.json",
        "--out-rev", tmp_path / "R.json", "--priors-fwd", tmp_path / "priors.fa", code=1,
    )  # fmt: skip
    assert "different numbers of reads" in res.stderr
    with open(tmp_path / "bad.pkl", "wb") as fh:
        pickle.dump([1, 2], fh)
    res = run(
        "infer-asv", "--in-fwd", d / "filt" / "sam1_F.fq.gz",
        "--in-rev", d / "filt" / "sam1_R.fq.gz",
        "--errors", tmp_path / "bad.pkl", "--out-fwd", tmp_path / "F.json",
        "--out-rev", tmp_path / "R.json", code=2,
    )  # fmt: skip
    assert "learn-errors" in res.stderr


# --------------------------------------------------------------------- merge errors ---
def test_merge_pairing_errors(tmp_path, pipeline):
    asv = pipeline["dir"] / "asv"
    common = ["--out-seqtab", tmp_path / "s.tsv", "--out-seqtab-nochim", tmp_path / "n.tsv"]
    res = run("merge", "--in-fwd", asv / "sam1_F.dada.json.gz", "--in-rev",
              asv / "sam2_R.dada.json.gz", *common, code=2)  # fmt: skip
    assert "sam1, sam2" in res.stderr
    res = run("merge", "--in-fwd", asv / "sam1_R.dada.json.gz", "--in-rev",
              asv / "sam1_R.dada.json.gz", *common, code=2)  # fmt: skip
    assert "expected 'fwd'" in res.stderr
    res = run("merge", "--in-fwd", asv / "sam1_F.dada.json.gz", asv / "sam1_F.dada.json.gz",
              "--in-rev", asv / "sam1_R.dada.json.gz", *common, code=2)  # fmt: skip
    assert "more than once" in res.stderr


@pytest.mark.parametrize("method", ["pooled", "per-sample"])
def test_merge_options(tmp_path, pipeline, method):
    asv = pipeline["dir"] / "asv"
    res = run(
        "merge", "--in-fwd", *[asv / f"{n}_F.dada.json.gz" for n in NAMES],
        "--in-rev", *[asv / f"{n}_R.dada.json.gz" for n in NAMES],
        "--out-seqtab", tmp_path / "s.tsv.gz", "--out-seqtab-nochim", tmp_path / "n.tsv",
        "--out-mergers", tmp_path / "m.tsv", "--return-rejects", "--order-by", "none",
        "--chimera-method", method, "--min-overlap", 300, "--stats", tmp_path / "st.tsv",
        code=EXIT_NO_DATA,
    )  # fmt: skip
    assert all(r["merged"] == "0" for r in tsv(res.stdout))
    m = pd.read_csv(tmp_path / "m.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert set(m["accept"]) == {"false"} and set(m["asv_id"]) == {""}
    res = run(
        "merge", "--in-fwd", *[asv / f"{n}_F.dada.json.gz" for n in NAMES],
        "--in-rev", *[asv / f"{n}_R.dada.json.gz" for n in NAMES],
        "--out-seqtab", tmp_path / "s.tsv.gz", "--out-seqtab-nochim", tmp_path / "n.tsv",
        "--just-concatenate", "--order-by", "nsamples", "--chimera-method", method,
    )  # fmt: skip
    seqtab = dada2.read_seqtab(tmp_path / "s.tsv.gz")
    assert all("NNNNNNNNNN" in s for s in seqtab.columns)


# ------------------------------------------------------------------ learn-errors etc ---
def test_learn_errors_options(tmp_path, pipeline):
    filt = pipeline["dir"] / "filt"
    f = [filt / f"{n}_F.fq.gz" for n in NAMES]
    r = [filt / f"{n}_R.fq.gz" for n in NAMES]
    res = run("learn-errors", "--in-fwd", *f, "--in-rev", r[0], "--out", tmp_path / "e.pkl",
              code=2)  # fmt: skip
    assert "as many files" in res.stderr
    base = ["learn-errors", "--in-fwd", *f, "--in-rev", *r, "--out", tmp_path / "e.pkl"]
    assert "requires --quality-bins" in run(*base, "--error-model", "binned", code=2).stderr
    assert "only used with" in run(*base, "--quality-bins", 2, 40, code=2).stderr
    assert "unsupported figure" in run(*base, "--plot", tmp_path / "x.abc", code=2).stderr
    res = run(*base, "--error-model", "noqual", "--max-consist", 1, "--randomize", "--seed", 1,
              "--n-bases", 1, "--band-size", 32, "--no-use-kmers", "--stats", tmp_path / "s.tsv",
              "--plot", tmp_path / "e.svg")  # fmt: skip
    rows = tsv(res.stdout)
    assert all(r["n_samples"] == "1" and r["converged"] == "false" for r in rows)
    assert "did not converge" in res.stderr
    expected = dada2.learn_errors(f, error_model=dada2.NoQualErrorModel(), max_consist=1,
                                  randomize=True, random_state=1, n_bases=1,
                                  params={"band_size": 32, "use_kmers": False})  # fmt: skip
    with open(tmp_path / "e.pkl", "rb") as fh:
        np.testing.assert_array_equal(pickle.load(fh)["err_f"], expected.err_)


def test_assign_taxonomy_options(tmp_path, pipeline):
    d = pipeline["dir"]
    res = run("assign-taxonomy", "--in", d / "nochim.tsv", "--ref", d / "model.pkl",
              "--out-taxa", tmp_path / "t.tsv", "--tax-levels", "K", "P", "C",
              "--min-boot", 0, "--try-rc", "--seed", 3, "--species-ref", SPECIES,
              "--allow-multiple", "--threads", 2)  # fmt: skip
    cols = (tmp_path / "t.tsv").read_text().split("\n")[0].split("\t")
    assert cols[:5] == ["asv_id", "sequence", "K", "P", "C"] and cols[-1] == "Species"
    assert tsv(res.stdout)[0]["rank"] == "K"
    run("assign-taxonomy", "--in", d / "nochim.tsv", "--ref", TRAIN, "--out-taxa",
        tmp_path / "t.tsv", "--species-ref", SPECIES, "--allow-multiple", 2)  # fmt: skip
    res = run("assign-taxonomy", "--in", d / "nochim.tsv", "--ref", TRAIN, "--out-taxa",
              tmp_path / "t.tsv", "--allow-multiple", 0, code=2)  # fmt: skip
    assert "positive integer" in res.stderr
    run("assign-taxonomy", "--in", d / "nochim.tsv", "--ref", TRAIN, "--out-taxa",
        tmp_path / "t.tsv", "--allow-multiple", "x", code=2)  # fmt: skip
    with open(tmp_path / "bad.pkl", "wb") as fh:
        pickle.dump({"a": 1}, fh)
    res = run("assign-taxonomy", "--in", d / "nochim.tsv", "--ref", tmp_path / "bad.pkl",
              "--out-taxa", tmp_path / "t.tsv", code=2)  # fmt: skip
    assert "TaxonomyClassifier" in res.stderr
    (tmp_path / "x.tsv").write_text("a\tb\n1\t2\n")
    res = run("assign-taxonomy", "--in", tmp_path / "x.tsv", "--ref", TRAIN,
              "--out-taxa", tmp_path / "t.tsv", code=2)  # fmt: skip
    assert "asv_id" in res.stderr


# ------------------------------------------------------------------------- formats ---
def test_dada_json_roundtrip_plain_and_errors(tmp_path):
    drp = dada2.derep_fastq(example_path("sam1F.fastq.gz"))
    res = dada2.DADA(err=err_balanced_f()).fit(drp).results_[0]
    dada2.write_dada_json(tmp_path / "x.json", res, drp, sample="s", mate="fwd")
    s = dada2.read_dada_json(tmp_path / "x.json")
    assert s.source is None
    pd.testing.assert_series_equal(s.result.denoised, res.denoised)
    assert s.result.quality.shape == (len(res.clusters), 0)
    # merge_pairs accepts the restored objects
    m = dada2.merge_pairs(s.result, s.derep, s.result, s.derep, min_overlap=0, max_mismatch=500)
    assert isinstance(m, pd.DataFrame)
    with pytest.raises(ValueError, match="mate"):
        dada2.write_dada_json(tmp_path / "y.json", res, drp, sample="s", mate="both")
    other = dada2.derep_fastq(example_path("sam2F.fastq.gz"))
    with pytest.raises(ValueError, match="does not correspond"):
        dada2.write_dada_json(tmp_path / "y.json", res, other, sample="s", mate="fwd")
    (tmp_path / "bad.json").write_text('{"format": "other"}')
    with pytest.raises(ValueError, match="format: 'other'"):
        dada2.read_dada_json(tmp_path / "bad.json")
    (tmp_path / "bad2.json").write_text("not json")
    with pytest.raises(ValueError, match=r"not a pydada2\.dada/1"):
        dada2.read_dada_json(tmp_path / "bad2.json")


def test_dada_json_is_reproducible(tmp_path):
    res, drp = _empty_dada(dada2.DadaParams())
    for name in ("a.json.gz", "b.json.gz"):
        dada2.write_dada_json(tmp_path / name, res, drp, sample="s", mate="rev")
    assert (tmp_path / "a.json.gz").read_bytes() == (tmp_path / "b.json.gz").read_bytes()


def test_seqtab_roundtrip(tmp_path):
    tab = pd.DataFrame({"ACGT": [1, 0], "GGCC": [5, 7]}, index=["s1", "s2"])
    dada2.write_seqtab(tab, tmp_path / "t.tsv.gz", asv_ids=["A", "B"])
    back, ids = dada2.read_seqtab(tmp_path / "t.tsv.gz", return_ids=True)
    pd.testing.assert_frame_equal(back, tab, check_column_type=False, check_index_type=False)
    assert ids.to_dict() == {"ACGT": "A", "GGCC": "B"}
    with pytest.raises(ValueError, match="one id per sequence"):
        dada2.write_seqtab(tab, tmp_path / "t.tsv", asv_ids=["A"])
    (tmp_path / "x.tsv").write_text("sequence\tasv_id\n")
    with pytest.raises(ValueError, match="asv_id"):
        dada2.read_seqtab(tmp_path / "x.tsv")
