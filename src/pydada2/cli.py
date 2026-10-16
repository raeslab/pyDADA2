"""Command-line interface: the DADA2 paired-end workflow as separate steps.

Each step is a subcommand of the ``pydada2`` console script (also ``python -m pydada2``), meant
to be orchestrated and parallelised by a workflow manager such as Nextflow::

    pydada2 trim-and-filter  ...   # per sample:  FASTQ pair -> filtered FASTQ pair
    pydada2 learn-errors     ...   # per batch:   filtered FASTQs -> error rates (pickle)
    pydada2 infer-asv        ...   # per sample:  filtered FASTQ pair -> 2 x .dada.json.gz
    pydada2 merge            ...   # per batch:   .dada.json.gz files -> sequence tables
    pydada2 assign-taxonomy  ...   # once:        sequence table -> taxonomy table

Conventions: options only (no positional arguments); per-mate parameters take ``FWD [REV]``;
file lists take several paths after one flag (``--in-fwd filt/*_F.fq.gz``); outputs are gzipped
when their name ends in ``.gz``; stdout carries only a TSV of read counts (also written to
``--stats``), logs go to stderr. Exit codes: 0 ok, 1 runtime error, 2 usage error, 3 no data.
"""

from __future__ import annotations

import logging
import math
import pickle
import traceback
import warnings
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, TextIO

import click
import numpy as np
import pandas as pd

from .chimeras import BimeraRemover
from .denoise import DADA, DadaParams, DadaResult
from .derep import Dereplication, derep_fastq
from .error_models import (
    TRANSITIONS,
    BinnedQualErrorModel,
    ErrorLearner,
    ErrorModel,
    LoessErrorModel,
    NoQualErrorModel,
    PacBioErrorModel,
    learn_errors,
)
from .filtering import filter_and_trim
from .formats import (
    DadaSample,
    _empty_dada,
    _open_write,
    read_dada_json,
    write_dada_json,
    write_seqtab,
)
from .io import QualityType, read_fastq_chunks, write_fasta
from .merge import merge_pairs
from .seqtab import make_sequence_table
from .taxonomy import DEFAULT_TAX_LEVELS, TaxonomyClassifier, add_species

if TYPE_CHECKING:  # pragma: no cover
    from click.parser import _OptionParser, _ParsingState

__all__ = ["main"]

logger = logging.getLogger("pydada2")

#: Exit code of a step that found no data (no read passed, no read merged, no ASV); its stats
#: are still written.
EXIT_NO_DATA = 3


# --------------------------------------------------------------------------- option types ---
def _is_option(token: str) -> bool:
    return token.startswith("-")


def _is_number(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


class GreedyOption(click.Option):
    """An option that takes one or more values after a single flag.

    ``--in-fwd a.fq b.fq c.fq`` gives ``("a.fq", "b.fq", "c.fq")``, so shell globs work. The
    values end at the first token that ``accept`` rejects (by default: anything starting with
    ``-``) or after ``max_values`` values. Repeating the flag adds more values. click has no
    variadic options, so the parser's handler of this option is wrapped to consume the
    following tokens.

    Parameters
    ----------
    *args, **kwargs
        As for :class:`click.Option` (``multiple`` is always True).
    max_values : int, optional
        Maximum number of values per occurrence of the flag.
    accept : callable, optional
        ``accept(token) -> bool``: whether a following token is another value.
    """

    def __init__(
        self,
        *args: Any,
        max_values: int | None = None,
        accept: Callable[[str], bool] | None = None,
        **kwargs: Any,
    ) -> None:
        kwargs["multiple"] = True
        super().__init__(*args, **kwargs)
        self.max_values = max_values
        self.accept = accept if accept is not None else (lambda tok: not _is_option(tok))

    def add_to_parser(self, parser: _OptionParser, ctx: click.Context) -> None:
        """Register the option, then make its parser handler consume the following values."""
        super().add_to_parser(parser, ctx)
        for name in self.opts:
            popt = parser._long_opt.get(name) or parser._short_opt.get(name)
            if popt is None:  # pragma: no cover - every opt is registered just above
                continue
            store = popt.process

            def process(value: Any, state: _ParsingState, store: Any = store) -> None:
                store(value, state)
                n = 1
                while (
                    state.rargs
                    and (self.max_values is None or n < self.max_values)
                    and self.accept(state.rargs[0])
                ):
                    store(state.rargs.pop(0), state)
                    n += 1

            popt.process = process  # type: ignore[method-assign]


class _IntOrInf(click.ParamType):  # type: ignore[type-arg]
    """An integer, or ``inf``/``-inf`` (thresholds that R compares with read properties)."""

    name = "integer|inf"

    def convert(self, value: Any, param: click.Parameter | None, ctx: click.Context | None) -> Any:
        if isinstance(value, int | float) and not isinstance(value, bool):
            return value
        text = str(value).strip().lower()
        if text in ("inf", "+inf", "-inf"):
            return -math.inf if text.startswith("-") else math.inf
        try:
            return int(text)
        except ValueError:
            self.fail(f"{value!r} is not an integer, 'inf' or '-inf'.", param, ctx)


class _AllowMultiple(click.ParamType):  # type: ignore[type-arg]
    """A positive integer, or ``all`` (the value of the bare flag)."""

    name = "N"

    def convert(self, value: Any, param: click.Parameter | None, ctx: click.Context | None) -> Any:
        if value is True or str(value).lower() == "all":
            return True
        try:
            n = int(value)
        except ValueError:
            self.fail(f"{value!r} is not a positive integer.", param, ctx)
        if n < 1:
            self.fail(f"{value!r} is not a positive integer.", param, ctx)
        return n


_IN_FILE = click.Path(exists=True, dir_okay=False, readable=True, path_type=Path)
_OUT_FILE = click.Path(dir_okay=False, writable=True, path_type=Path)
_QUALITY_TYPES = click.Choice(["auto", "FastqQuality", "SFastqQuality"])

Decorator = Callable[[Callable[..., Any]], Callable[..., Any]]


def _compose(*decorators: Decorator) -> Decorator:
    """Apply decorators so that options appear in ``--help`` in the order given."""

    def wrap(f: Callable[..., Any]) -> Callable[..., Any]:
        for d in reversed(decorators):
            f = d(f)
        return f

    return wrap


def _file_list(flag: str, help: str) -> Decorator:
    return click.option(
        flag, cls=GreedyOption, type=_IN_FILE, required=True, metavar="PATH...", help=help
    )


def _per_mate(flag: str, type_: Any, default: Any, help: str) -> Decorator:
    shown = "inf" if isinstance(default, float) and math.isinf(default) else str(default)
    return click.option(
        flag,
        cls=GreedyOption,
        type=type_,
        max_values=2,
        accept=_is_number,
        default=(default,),
        show_default=False,
        metavar="FWD [REV]",
        callback=_mate_values,
        help=f"{help} One value for both mates, or two values (forward, reverse).  "
        f"[default: {shown}]",
    )


def _mate_values(ctx: click.Context, param: click.Parameter, value: tuple[Any, ...]) -> Any:
    if len(value) == 1:
        return value[0]
    if len(value) == 2:
        return tuple(value)
    raise click.BadParameter("give one value (both mates) or two values (forward, reverse).")


_common_options = _compose(
    click.option(
        "--threads",
        type=click.IntRange(min=1),
        default=1,
        show_default=True,
        help="Number of worker threads.",
    ),
    click.option(
        "-v",
        "--verbose",
        count=True,
        help="Log progress to stderr (repeat for more detail: -vv also shows tracebacks).",
    ),
)

_stats_option = click.option(
    "--stats",
    type=_OUT_FILE,
    default=None,
    help="Also write the read-tracking TSV printed on stdout to this file.",
)

_sample_name_option = click.option(
    "--sample-name",
    metavar="NAME",
    default=None,
    help="Sample name. Default: the --in-fwd file name without .gz, .fastq/.fq/.dada.json and "
    "a trailing _R1/_1/_F (sam1_R1.fastq.gz -> sam1).",
)


def _quality_type_option(help_suffix: str = "") -> Decorator:
    return click.option(
        "--quality-type",
        type=_QUALITY_TYPES,
        default="auto",
        show_default=True,
        help="Quality encoding of the FASTQ files: auto-detected, FastqQuality (Phred+33) or "
        "SFastqQuality (Solexa/Phred+64)." + help_suffix,
    )


# DADA algorithm options (R: setDadaOpt), one CLI option per DadaParams field.
_DADA_HELP: dict[str, str] = {
    "omega_a": "P-value threshold of the abundance test for calling a unique sequence "
    "over-abundant, which makes it the centre of a new partition (ASV). Lower is more "
    "conservative.",
    "omega_p": "P-value threshold for unique sequences with prior evidence of existence "
    "(--priors-fwd/--priors-rev).",
    "omega_c": "Error-correction threshold: a unique sequence is replaced by its partition "
    "centre if its p-value is >= this, otherwise it is left uncorrected (dropped from the "
    "denoised counts). 0 corrects everything.",
    "detect_singletons": "Allow sequences seen in a single read to form new ASVs (more "
    "sensitive to rare variants, more false positives).",
    "use_kmers": "Screen pairs of sequences by 5-mer distance before aligning them.",
    "kdist_cutoff": "5-mer distance above which a pair of sequences is not aligned (with "
    "--use-kmers).",
    "match": "Needleman-Wunsch alignment score of a match.",
    "mismatch": "Needleman-Wunsch alignment score of a mismatch.",
    "gap_penalty": "Needleman-Wunsch alignment score of a gap.",
    "homopolymer_gap_penalty": "Alignment score of gaps in homopolymers (>= 3 repeated "
    "bases). Default: the same as --gap-penalty.",
    "band_size": "Band width of the banded alignment (net insertions of one sequence relative "
    "to the other); negative: unbanded. Increase for indel-rich technologies.",
    "max_clust": "Maximum number of partitions (ASVs) per sample; 0: no limit.",
    "min_fold": "Minimum fold over-abundance (observed / expected) for a new partition.",
    "min_hamming": "Minimum Hamming distance to the partition centre for a new partition.",
    "min_abundance": "Minimum abundance of a unique sequence to form a new partition.",
    "use_quals": "Use quality scores in the error model (off: quality-independent rates).",
    "gapless": "Assume a gapless alignment when the ordered k-mer identity of two sequences "
    "equals their unordered overlap (with --use-kmers).",
    "greedy": "Lock sequences with fewer reads than expected from resequencing the partition "
    "centre to that partition (~30% faster, almost no effect on the output).",
    "vectorized_alignment": "Use the vectorized (SIMD) aligner.",
    "sse": "SIMD level of the k-mer routines: 0 (off), 1 or 2.",
}
# Not exposed: max_consist (learn-errors has its own option; infer-asv does not iterate) and
# pseudo_prevalence/pseudo_abundance (pseudo-pooling is not a per-sample step).


def _dada_options(exclude: Sequence[str] = ()) -> Decorator:
    defaults = DadaParams()
    decorators: list[Decorator] = []
    for name, help_ in _DADA_HELP.items():
        if name in exclude:
            continue
        flag = "--" + name.replace("_", "-")
        default = getattr(defaults, name)
        if isinstance(default, bool):
            decorators.append(
                click.option(
                    f"{flag}/--no-{name.replace('_', '-')}",
                    name,
                    default=default,
                    show_default=True,
                    help=help_,
                )
            )
        elif name == "homopolymer_gap_penalty":
            decorators.append(click.option(flag, name, type=int, default=None, help=help_))
        else:
            decorators.append(
                click.option(
                    flag, name, type=type(default), default=default, show_default=True, help=help_
                )
            )
    return _compose(*decorators)


def _pop_dada_params(kwargs: dict[str, Any]) -> dict[str, Any]:
    return {k: kwargs.pop(k) for k in list(kwargs) if k in _DADA_HELP}


# ---------------------------------------------------------------------------- helpers ---
class _StderrHandler(logging.Handler):
    """Logging handler writing to click's (current) stderr."""

    def emit(self, record: logging.LogRecord) -> None:
        try:
            click.echo(self.format(record), err=True)
        except Exception:  # pragma: no cover - logging must never raise
            self.handleError(record)


def _show_warning(
    message: Warning | str,
    category: type[Warning],
    filename: str,
    lineno: int,
    file: TextIO | None = None,
    line: str | None = None,
) -> None:
    click.echo(f"Warning: {message}", err=True)


@contextmanager
def _session(verbose: int) -> Iterator[None]:
    """Run a subcommand: logs and warnings to stderr; errors -> message and exit code 1."""
    handler = _StderrHandler()
    handler.setFormatter(logging.Formatter("%(message)s"))
    old_level = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose >= 2 else logging.INFO if verbose else logging.WARNING)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("default")
            warnings.showwarning = _show_warning
            try:
                yield
            except (click.ClickException, click.exceptions.Exit, click.Abort):
                raise
            except Exception as err:
                if verbose >= 2:
                    click.echo(traceback.format_exc(), err=True, nl=False)
                raise click.ClickException(str(err) or type(err).__name__) from err
    finally:
        logger.removeHandler(handler)
        logger.setLevel(old_level)


def sample_name_from_path(path: Path | str) -> str:
    """Derive a sample name from a file name.

    Strips ``.gz``, then ``.fastq``, ``.fq`` or ``.dada.json``, then a trailing ``_R1``,
    ``_1`` or ``_F``.

    Parameters
    ----------
    path : str or pathlib.Path
        File path.

    Returns
    -------
    str
        The sample name.

    Examples
    --------
    >>> sample_name_from_path("reads/sam1_R1.fastq.gz")
    'sam1'
    >>> sample_name_from_path("asv/sam2_F.dada.json.gz")
    'sam2'
    """
    name = Path(path).name
    name = name.removesuffix(".gz")
    for ext in (".dada.json", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name.removesuffix(ext)
            break
    for suffix in ("_R1", "_1", "_F"):
        if name.endswith(suffix) and len(name) > len(suffix):
            name = name.removesuffix(suffix)
            break
    return name


def _prepare_outputs(inputs: Sequence[Path], outputs: Sequence[Path | None]) -> None:
    """Check that outputs are distinct from each other and from the inputs; create dirs."""
    outs = [p for p in outputs if p is not None]
    resolved = [p.resolve() for p in outs]
    if len(set(resolved)) != len(resolved):
        raise click.UsageError("All output files must be distinct.")
    ins = {p.resolve() for p in inputs}
    for p, r in zip(outs, resolved, strict=True):
        if r in ins:
            raise click.UsageError(f"Output file {p} is also an input file.")
    for p in outs:
        p.parent.mkdir(parents=True, exist_ok=True)


def _fmt(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)) or value is pd.NA:
        return ""
    if isinstance(value, bool | np.bool_):
        return "true" if value else "false"
    return str(value)


def _emit_stats(rows: Sequence[dict[str, Any]], columns: Sequence[str], path: Path | None) -> None:
    """Print a TSV (header + rows) on stdout and, optionally, write it to ``path``."""
    lines = ["\t".join(columns)]
    lines += ["\t".join(_fmt(r[c]) for c in columns) for r in rows]
    text = "\n".join(lines) + "\n"
    click.echo(text, nl=False)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with _open_write(path) as fh:
            fh.write(text.encode())


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    out = df.copy()
    for c in out.columns:
        if out[c].dtype == bool:
            out[c] = out[c].map(_fmt)
        elif out[c].dtype == np.float64:  # shortest repr that reads back exactly
            out[c] = [_fmt(v) if math.isnan(v) else repr(float(v)) for v in out[c]]
    with _open_write(path) as fh:
        out.to_csv(fh, sep="\t", index=False, na_rep="")


def _no_data() -> None:
    raise click.exceptions.Exit(EXIT_NO_DATA)


def _has_reads(path: Path) -> bool:
    return next(iter(read_fastq_chunks(path, n=1)), None) is not None


# ------------------------------------------------------------------------------- group ---
@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 100},
)
@click.version_option(package_name="pydada2", prog_name="pydada2")
def cli() -> None:
    """DADA2 amplicon workflow (paired-end reads) as separate steps.

    \b
    Typical order:
      trim-and-filter   per sample   raw FASTQ pair -> filtered FASTQ pair
      learn-errors      per batch    filtered FASTQs -> error rates (errors.pkl)
      infer-asv         per sample   filtered FASTQ pair -> denoised mates (.dada.json.gz)
      merge             per batch    denoised mates -> sequence tables (+ chimera removal)
      assign-taxonomy   once         sequence table -> taxonomy table

    stdout carries only a TSV of read counts (header line + rows); logs and warnings go to
    stderr. Exit codes: 0 ok, 1 runtime error, 2 usage error, 3 no data (stats still written).
    Please cite the DADA2 paper (Callahan et al. 2016, Nature Methods) when using pydada2.
    """  # noqa: D301 - click: \b keeps line breaks; \\ shows a backslash


# --------------------------------------------------------------------- trim-and-filter ---
@cli.command("trim-and-filter")
@click.option(
    "--in-fwd", type=_IN_FILE, required=True, help="Forward-read FASTQ file (optionally gzipped)."
)
@click.option(
    "--in-rev",
    type=_IN_FILE,
    required=True,
    help="Matching reverse-read FASTQ file (same reads, same order).",
)
@click.option(
    "--out-fwd",
    type=_OUT_FILE,
    required=True,
    help="Filtered forward reads (FASTQ; gzipped if the name ends in .gz). Not "
    "written if no read passes.",
)
@click.option(
    "--out-rev",
    type=_OUT_FILE,
    required=True,
    help="Filtered reverse reads (compressed like --out-fwd).",
)
@_per_mate("--trunc-q", int, 2, "Truncate reads at the first base with a quality score <= this.")
@_per_mate(
    "--trunc-len",
    int,
    0,
    "Truncate reads after this many bases; shorter reads are discarded (0: off).",
)
@_per_mate("--trim-left", int, 0, "Number of bases removed from the start of each read.")
@_per_mate("--trim-right", int, 0, "Number of bases removed from the end of each read.")
@_per_mate(
    "--max-len", float, math.inf, "Discard reads longer than this, before trimming and truncation."
)
@_per_mate(
    "--min-len", _IntOrInf(), 20, "Discard reads shorter than this, after trimming and truncation."
)
@_per_mate(
    "--max-n", _IntOrInf(), 0, "Maximum number of ambiguous (non-ACGT) bases; DADA2 requires 0."
)
@_per_mate(
    "--min-q",
    _IntOrInf(),
    0,
    "Discard reads with any quality score <= this, after truncation (only applied when "
    "greater than --trunc-q, as in R).",
)
@_per_mate(
    "--max-ee",
    float,
    math.inf,
    "Discard reads with more expected errors, sum(10^(-Q/10)), than this.",
)
@_per_mate(
    "--rm-lowcomplex",
    float,
    0,
    "Discard reads with a sequence complexity (effective number of 2-mers) below this (0: off).",
)
@click.option(
    "--rm-phix/--keep-phix",
    default=True,
    show_default=True,
    help="Remove read pairs matching the phiX genome (Illumina spike-in).",
)
@click.option(
    "--orient-fwd",
    metavar="SEQ",
    default=None,
    help="Keep only pairs whose forward or reverse read starts with SEQ (e.g. the "
    "forward primer); mates are swapped where the reverse read matches.",
)
@click.option(
    "--match-ids/--no-match-ids",
    default=False,
    show_default=True,
    help="Keep only pairs whose read ids match in both files (the files must be in "
    "the same relative order).",
)
@click.option(
    "--id-sep",
    metavar="REGEX",
    default=r"\s",
    show_default=True,
    help="Regular expression separating the fields of the read ids (--match-ids).",
)
@click.option(
    "--id-field",
    type=click.IntRange(min=0),
    default=None,
    help="0-based id field compared by --match-ids. Default: auto-detected (Illumina "
    "CASAVA formats).",
)
@_quality_type_option()
@click.option(
    "--chunk-size",
    type=click.IntRange(min=1),
    default=100_000,
    show_default=True,
    help="Reads processed per chunk (memory use; affects the output order with "
    "--orient-fwd, as in R).",
)
@_sample_name_option
@_stats_option
@_common_options
def trim_and_filter(
    in_fwd: Path,
    in_rev: Path,
    out_fwd: Path,
    out_rev: Path,
    orient_fwd: str | None,
    match_ids: bool,
    id_sep: str,
    id_field: int | None,
    quality_type: QualityType,
    chunk_size: int,
    sample_name: str | None,
    stats: Path | None,
    threads: int,
    verbose: int,
    **mate_params: Any,
) -> None:
    """Filter and trim one pair of FASTQ files (R: filterAndTrim).

    A read pair is kept only if both mates pass. Prints one TSV row (sample, reads_in,
    reads_out) and exits with code 3 when no read passes (no FASTQ is written then).

    \b
    Example:
      pydada2 trim-and-filter --in-fwd sam1_R1.fastq.gz --in-rev sam1_R2.fastq.gz \\
        --out-fwd filtered/sam1_F.fq.gz --out-rev filtered/sam1_R.fq.gz \\
        --trunc-len 240 160 --max-ee 2 --trunc-q 2 --stats sam1.filter_stats.tsv
    """  # noqa: D301 - click: \b keeps line breaks; \\ shows a backslash
    with _session(verbose):
        _prepare_outputs([in_fwd, in_rev], [out_fwd, out_rev])
        compress = out_fwd.name.endswith(".gz")
        if out_rev.name.endswith(".gz") != compress:
            raise click.UsageError("--out-fwd and --out-rev must both end in .gz, or neither.")
        name = sample_name or sample_name_from_path(in_fwd)
        with warnings.catch_warnings():  # reported below, with the sample name
            warnings.filterwarnings("ignore", message="No reads passed the filter")
            df = filter_and_trim(
                in_fwd,
                out_fwd,
                in_rev,
                out_rev,
                compress=compress,
                orient_fwd=orient_fwd,
                match_ids=match_ids,
                id_sep=id_sep,
                id_field=id_field,
                n_jobs=threads,
                n=chunk_size,
                quality_type=quality_type,
                verbose=verbose,
                **mate_params,
            )
        reads_in, reads_out = int(df["reads_in"].iloc[0]), int(df["reads_out"].iloc[0])
        row = {"sample": name, "reads_in": reads_in, "reads_out": reads_out}
        _emit_stats([row], ["sample", "reads_in", "reads_out"], stats)
        if reads_out == 0:
            logger.warning("%s: no reads passed the filter; no output written.", name)
            _no_data()


# ------------------------------------------------------------------------ learn-errors ---
def _error_model(name: str, quality_bins: Sequence[int]) -> ErrorModel:
    if quality_bins and name != "binned":
        raise click.UsageError("--quality-bins is only used with --error-model binned.")
    if name == "binned":
        if not quality_bins:
            raise click.UsageError("--error-model binned requires --quality-bins.")
        return BinnedQualErrorModel(tuple(quality_bins))
    models: dict[str, Callable[[], ErrorModel]] = {
        "loess": LoessErrorModel,
        "noqual": NoQualErrorModel,
        "pacbio": PacBioErrorModel,
    }
    return models[name]()


def _check_plot_path(path: Path) -> None:
    try:
        from matplotlib.backend_bases import FigureCanvasBase
    except ImportError:  # pragma: no cover - depends on the environment
        raise click.UsageError("--plot requires matplotlib (pip install pydada2[plot]).") from None
    ext = path.suffix.lstrip(".").lower()
    if ext not in FigureCanvasBase.get_supported_filetypes():
        raise click.BadParameter(
            f"unsupported figure format {ext!r} (use e.g. .png, .pdf or .svg).",
            param_hint="--plot",
        )


def _plot_errors_pair(learners: dict[str, ErrorLearner], path: Path) -> None:
    """Draw R's plotErrors panels for both mates side by side and save the figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .plotting import error_plot_data

    fig = plt.figure(figsize=(12 * len(learners), 10), layout="constrained")
    subfigs = np.atleast_1d(fig.subfigures(1, len(learners)))  # type: ignore[type-var]
    for sub, (mate, learner) in zip(subfigs, learners.items(), strict=True):  # type: ignore[arg-type]
        df = error_plot_data(learner)
        sub.suptitle(f"{mate} reads", fontsize="x-large")
        axes = sub.subplots(4, 4, squeeze=False)
        for i, a in enumerate("ACGT"):
            for j, b in enumerate("ACGT"):  # as plot_errors(): observed points, estimated line
                ax = axes[i][j]
                part = df[df["transition"] == f"{a}2{b}"]
                if part["observed"].notna().any():
                    cnt = part["count"].fillna(0).to_numpy()
                    size = 16 * np.sqrt(cnt / cnt.max()) if cnt.max() > 0 else 1
                    ax.scatter(part["qual"], part["observed"], s=size, color="0.4")
                if part["estimated"].notna().any():
                    ax.plot(part["qual"], part["estimated"], color="black")
                ax.set_yscale("log")
                ax.set_title(f"{a}2{b}")
                if i == 3:
                    ax.set_xlabel("Consensus quality score")
                if j == 0:
                    ax.set_ylabel("Error frequency (log10)")
    fig.savefig(path)
    plt.close(fig)


@cli.command("learn-errors")
@_file_list(
    "--in-fwd", "Filtered forward FASTQ files (several paths after one flag; shell globs work)."
)
@_file_list(
    "--in-rev", "Filtered reverse FASTQ files (as many as --in-fwd; order does not matter)."
)
@click.option(
    "--out",
    type=_OUT_FILE,
    required=True,
    help="Error rates: a pickle of {'err_f': ndarray, 'err_r': ndarray}, the 16 x Q "
    "error matrices (rows A2A, A2C, ..., T2T; columns quality 0..Q-1).",
)
@click.option(
    "--out-tsv",
    type=_OUT_FILE,
    default=None,
    help="Also write the error rates as a long TSV: mate, transition, qual, rate.",
)
@click.option(
    "--plot",
    type=_OUT_FILE,
    default=None,
    help="Plot observed and fitted error rates of both mates (R: plotErrors); the "
    "format follows the extension (png, pdf, svg, ...). Needs matplotlib.",
)
@click.option(
    "--n-bases",
    type=click.FloatRange(min=0),
    default=1e8,
    show_default=True,
    help="Minimum number of bases to learn from, per mate: files are read in order "
    "until this many bases have been read.",
)
@click.option(
    "--error-model",
    type=click.Choice(["loess", "noqual", "binned", "pacbio"]),
    default="loess",
    show_default=True,
    help="How error rates depend on quality: loess fit (default), constant (noqual), "
    "interpolated between binned quality scores (binned, e.g. NovaSeq), or PacBio "
    "CCS.",
)
@click.option(
    "--quality-bins",
    cls=GreedyOption,
    type=int,
    metavar="Q...",
    default=(),
    help="The quality scores present in binned data, e.g. 2 12 23 37 (with --error-model binned).",
)
@click.option(
    "--randomize/--no-randomize",
    default=False,
    show_default=True,
    help="Pick the files in random order instead of the given order.",
)
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for --randomize (default: not reproducible).",
)
@click.option(
    "--max-consist",
    type=click.IntRange(min=1),
    default=10,
    show_default=True,
    help="Maximum number of self-consistency iterations (denoise, re-estimate the "
    "error rates, repeat).",
)
@click.option(
    "--omega-c",
    "omega_c",
    type=float,
    default=0.0,
    show_default=True,
    help="Error-correction threshold while learning (see infer-asv); 0 corrects all "
    "reads, as R's learnErrors.",
)
@_dada_options(exclude=("omega_c",))
@_quality_type_option()
@_stats_option
@_common_options
def learn_errors_cmd(
    in_fwd: tuple[Path, ...],
    in_rev: tuple[Path, ...],
    out: Path,
    out_tsv: Path | None,
    plot: Path | None,
    n_bases: float,
    error_model: str,
    quality_bins: tuple[int, ...],
    randomize: bool,
    seed: int | None,
    max_consist: int,
    omega_c: float,
    quality_type: QualityType,
    stats: Path | None,
    threads: int,
    verbose: int,
    **dada_opts: Any,
) -> None:
    """Learn the error rates of the forward and of the reverse reads (R: learnErrors).

    Prints one TSV row per mate: mate, n_samples, n_reads, n_bases (used for learning),
    n_iter (self-consistency rounds) and converged. Not converging is a warning, not an error.

    \b
    Example:
      pydada2 learn-errors --in-fwd filtered/*_F.fq.gz --in-rev filtered/*_R.fq.gz \\
        --out errors.pkl --out-tsv errors.tsv.gz --plot errors.pdf --threads 8
    """  # noqa: D301 - click: \b keeps line breaks; \\ shows a backslash
    with _session(verbose):
        if len(in_fwd) != len(in_rev):
            raise click.UsageError(
                f"--in-fwd and --in-rev must list as many files ({len(in_fwd)} vs {len(in_rev)})."
            )
        _prepare_outputs([*in_fwd, *in_rev], [out, out_tsv, plot, stats])
        model = _error_model(error_model, quality_bins)
        if plot is not None:
            _check_plot_path(plot)
        params = _pop_dada_params(dada_opts)
        learners: dict[str, ErrorLearner] = {}
        for mate, files in (("fwd", in_fwd), ("rev", in_rev)):
            logger.info("Learning the error rates of the %s reads.", mate)
            learners[mate] = learn_errors(
                list(files),
                n_bases=n_bases,
                error_model=model,
                randomize=randomize,
                random_state=seed,
                max_consist=max_consist,
                omega_c=omega_c,
                params=params,
                quality_type=quality_type,
                n_jobs=threads,
                verbose=verbose,
            )
        err: dict[str, Any] = {"err_f": learners["fwd"].err_, "err_r": learners["rev"].err_}
        with _open_write(out) as fh:
            pickle.dump(err, fh, protocol=pickle.HIGHEST_PROTOCOL)
        if out_tsv is not None:
            long = [
                {"mate": mate, "transition": t, "qual": q, "rate": float(m[i, q])}
                for mate, m in (("fwd", err["err_f"]), ("rev", err["err_r"]))
                for i, t in enumerate(TRANSITIONS)
                for q in range(m.shape[1])
            ]
            _write_tsv(pd.DataFrame(long), out_tsv)
        if plot is not None:
            _plot_errors_pair(learners, plot)
        rows = []
        for mate, lrn in learners.items():
            if not lrn.converged_:
                logger.warning(
                    "The error rates of the %s reads did not converge in %d rounds (consider "
                    "a larger --max-consist).",
                    mate,
                    lrn.n_iter_,
                )
            rows.append(
                {
                    "mate": mate,
                    "n_samples": lrn.n_samples_,
                    "n_reads": lrn.n_reads_,
                    "n_bases": lrn.n_bases_,
                    "n_iter": lrn.n_iter_,
                    "converged": lrn.converged_,
                }
            )
        _emit_stats(rows, ["mate", "n_samples", "n_reads", "n_bases", "n_iter", "converged"], stats)


# --------------------------------------------------------------------------- infer-asv ---
def _load_errors(path: Path) -> tuple[Any, Any]:
    with open(path, "rb") as raw:
        gz = raw.read(2) == b"\x1f\x8b"
    if gz:
        import gzip

        with gzip.open(path, "rb") as fh:
            obj = pickle.load(fh)  # trusted: written by learn-errors
    else:
        with open(path, "rb") as fh:
            obj = pickle.load(fh)
    if not isinstance(obj, dict) or not {"err_f", "err_r"} <= set(obj):
        raise click.BadParameter(
            "expected a pickle of {'err_f': ..., 'err_r': ...} written by learn-errors.",
            param_hint="--errors",
        )
    return obj["err_f"], obj["err_r"]


@cli.command("infer-asv")
@click.option("--in-fwd", type=_IN_FILE, required=True, help="Filtered forward FASTQ file.")
@click.option(
    "--in-rev",
    type=_IN_FILE,
    required=True,
    help="Filtered reverse FASTQ file (same read pairs as --in-fwd).",
)
@click.option(
    "--errors",
    type=_IN_FILE,
    required=True,
    help="Error rates written by learn-errors (errors.pkl).",
)
@click.option(
    "--out-fwd",
    type=_OUT_FILE,
    required=True,
    help="Denoised forward reads (pydada2.dada/1 JSON, e.g. sam1_F.dada.json.gz).",
)
@click.option(
    "--out-rev", type=_OUT_FILE, required=True, help="Denoised reverse reads (pydada2.dada/1 JSON)."
)
@_sample_name_option
@click.option(
    "--priors-fwd",
    type=_IN_FILE,
    default=None,
    metavar="FASTA",
    help="Sequences with prior evidence of existence in the forward reads (tested "
    "against --omega-p).",
)
@click.option(
    "--priors-rev",
    type=_IN_FILE,
    default=None,
    metavar="FASTA",
    help="Prior sequences for the reverse reads.",
)
@click.option(
    "--chunk-size",
    type=click.IntRange(min=1),
    default=1_000_000,
    show_default=True,
    help="Reads per chunk while dereplicating (memory use; as in R, the order of "
    "equally abundant uniques depends on it).",
)
@_quality_type_option()
@_dada_options()
@_stats_option
@_common_options
def infer_asv(
    in_fwd: Path,
    in_rev: Path,
    errors: Path,
    out_fwd: Path,
    out_rev: Path,
    sample_name: str | None,
    priors_fwd: Path | None,
    priors_rev: Path | None,
    chunk_size: int,
    quality_type: QualityType,
    stats: Path | None,
    threads: int,
    verbose: int,
    **dada_opts: Any,
) -> None:
    """Dereplicate and denoise the forward and the reverse reads of one sample (R: dada).

    Does not merge the mates: each output keeps what merge needs (clusters, read maps, ...).
    Prints one TSV row: sample, reads, uniques_fwd, uniques_rev, denoised_fwd, denoised_rev.
    With zero input reads, both outputs are written with empty tables and the exit code is 3.

    \b
    Example:
      pydada2 infer-asv --in-fwd filtered/sam1_F.fq.gz --in-rev filtered/sam1_R.fq.gz \\
        --errors errors.pkl --out-fwd asv/sam1_F.dada.json.gz --out-rev asv/sam1_R.dada.json.gz
    """  # noqa: D301 - click: \b keeps line breaks; \\ shows a backslash
    with _session(verbose):
        _prepare_outputs([in_fwd, in_rev, errors], [out_fwd, out_rev, stats])
        name = sample_name or sample_name_from_path(in_fwd)
        err_f, err_r = _load_errors(errors)
        params = _pop_dada_params(dada_opts)
        columns = ["sample", "reads", "uniques_fwd", "uniques_rev", "denoised_fwd", "denoised_rev"]
        if not _has_reads(in_fwd) and not _has_reads(in_rev):
            empty_result, empty_derep = _empty_dada(DadaParams().replace(**params).resolved())
            for path, mate, src in ((out_fwd, "fwd", in_fwd), (out_rev, "rev", in_rev)):
                write_dada_json(
                    path, empty_result, empty_derep, sample=name, mate=mate, source=src.name
                )
            _emit_stats([dict.fromkeys(columns, 0) | {"sample": name}], columns, stats)
            logger.warning("%s: no input reads.", name)
            _no_data()

        results: dict[str, tuple[DadaResult, Dereplication]] = {}
        for mate, path, err, priors in (
            ("fwd", in_fwd, err_f, priors_fwd),
            ("rev", in_rev, err_r, priors_rev),
        ):
            drp = derep_fastq(
                path,
                n=chunk_size,
                quality_type=quality_type,
                verbose=verbose,
            )
            dd = DADA(
                err=err,
                priors=priors if priors is not None else (),
                params=params,
                n_jobs=threads,
                verbose=verbose,
            ).fit(drp)
            results[mate] = (dd.results_[0], drp)
        (res_f, drp_f), (res_r, drp_r) = results["fwd"], results["rev"]
        if len(drp_f.map) != len(drp_r.map):
            raise click.ClickException(
                f"{in_fwd.name} and {in_rev.name} have different numbers of reads "
                f"({len(drp_f.map)} vs {len(drp_r.map)}); they must hold the same read pairs."
            )
        write_dada_json(out_fwd, res_f, drp_f, sample=name, mate="fwd", source=in_fwd.name)
        write_dada_json(out_rev, res_r, drp_r, sample=name, mate="rev", source=in_rev.name)
        row = {
            "sample": name,
            "reads": len(drp_f.map),
            "uniques_fwd": len(drp_f),
            "uniques_rev": len(drp_r),
            "denoised_fwd": int(res_f.denoised.sum()),
            "denoised_rev": int(res_r.denoised.sum()),
        }
        _emit_stats([row], columns, stats)


# ------------------------------------------------------------------------------- merge ---
def _load_mates(paths: Sequence[Path], mate: str, flag: str) -> dict[str, DadaSample]:
    out: dict[str, DadaSample] = {}
    for p in paths:
        s = read_dada_json(p)
        if s.mate != mate:
            raise click.BadParameter(
                f"{p} holds {s.mate!r} reads, expected {mate!r}.", param_hint=flag
            )
        if s.sample in out:
            raise click.BadParameter(f"sample {s.sample!r} occurs more than once.", param_hint=flag)
        out[s.sample] = s
    return out


@cli.command("merge")
@_file_list("--in-fwd", "Denoised forward reads (.dada.json.gz files from infer-asv).")
@_file_list(
    "--in-rev",
    "Denoised reverse reads; paired with --in-fwd by the stored sample name (not by position).",
)
@click.option(
    "--out-seqtab",
    type=_OUT_FILE,
    required=True,
    help="Sequence table before chimera removal (TSV: asv_id, sequence, one column "
    "of counts per sample).",
)
@click.option(
    "--out-seqtab-nochim",
    type=_OUT_FILE,
    required=True,
    help="Sequence table after chimera removal (same ASV ids as --out-seqtab).",
)
@click.option(
    "--out-fasta",
    type=_OUT_FILE,
    default=None,
    help="Non-chimeric ASV sequences as FASTA (>ASV1 ...).",
)
@click.option(
    "--out-mergers",
    type=_OUT_FILE,
    default=None,
    help="Long TSV of every merged pair per sample (R: mergePairs output; forward/"
    "reverse are 0-based ASV indices of the denoised mates).",
)
@click.option(
    "--min-overlap",
    type=click.IntRange(min=0),
    default=12,
    show_default=True,
    help="Minimum overlap (in bases) of the forward and reverse reads.",
)
@click.option(
    "--max-mismatch",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    help="Maximum number of mismatches (and indels) allowed in the overlap.",
)
@click.option(
    "--return-rejects/--no-return-rejects",
    default=False,
    show_default=True,
    help="Also list the pairs that could not be merged in --out-mergers (accept = false).",
)
@click.option(
    "--just-concatenate/--no-just-concatenate",
    default=False,
    show_default=True,
    help="Concatenate the reads with 10 Ns in between instead of merging them (for "
    "non-overlapping reads).",
)
@click.option(
    "--trim-overhang/--no-trim-overhang",
    default=False,
    show_default=True,
    help="Trim the parts of the reads that extend past the start of their mate.",
)
@click.option(
    "--order-by",
    type=click.Choice(["abundance", "nsamples", "none"]),
    default="abundance",
    show_default=True,
    help="Order (and numbering) of the ASVs: by total abundance, by number of "
    "samples, or by first appearance.",
)
@click.option(
    "--chimera-method",
    type=click.Choice(["consensus", "pooled", "per-sample"]),
    default="consensus",
    show_default=True,
    help="consensus: flag chimeras per sample and remove those flagged in most "
    "samples; pooled: test the pooled counts; per-sample: remove per sample.",
)
@click.option(
    "--min-fold-parent-over-abundance",
    type=click.FloatRange(min=0),
    default=None,
    help="Parents must be this many times more abundant than the chimera. Default: "
    "1.5 for consensus, 2 otherwise (as in R).",
)
@click.option(
    "--min-parent-abundance",
    type=click.FloatRange(min=0),
    default=None,
    help="Minimum abundance of a parent. Default: 2 for consensus, 8 otherwise.",
)
@click.option(
    "--min-sample-fraction",
    type=click.FloatRange(0, 1),
    default=0.9,
    show_default=True,
    help="Fraction of the samples (containing the ASV) that must flag it as chimeric "
    "(consensus only).",
)
@click.option(
    "--ignore-n-negatives",
    type=click.IntRange(min=0),
    default=1,
    show_default=True,
    help="Samples that do not flag the ASV that are ignored in that fraction (consensus only).",
)
@click.option(
    "--allow-one-off/--no-allow-one-off",
    default=False,
    show_default=True,
    help="Also flag sequences one mismatch away from a two-parent chimera.",
)
@click.option(
    "--min-one-off-parent-distance",
    type=click.IntRange(min=0),
    default=4,
    show_default=True,
    help="Minimum distance of the parents from the sequence (with --allow-one-off).",
)
@click.option(
    "--max-shift",
    type=click.IntRange(min=0),
    default=16,
    show_default=True,
    help="Maximum shift when aligning the sequence to its potential parents.",
)
@_stats_option
@_common_options
def merge(
    in_fwd: tuple[Path, ...],
    in_rev: tuple[Path, ...],
    out_seqtab: Path,
    out_seqtab_nochim: Path,
    out_fasta: Path | None,
    out_mergers: Path | None,
    min_overlap: int,
    max_mismatch: int,
    return_rejects: bool,
    just_concatenate: bool,
    trim_overhang: bool,
    order_by: str,
    chimera_method: str,
    min_fold_parent_over_abundance: float | None,
    min_parent_abundance: float | None,
    min_sample_fraction: float,
    ignore_n_negatives: int,
    allow_one_off: bool,
    min_one_off_parent_distance: int,
    max_shift: int,
    stats: Path | None,
    threads: int,
    verbose: int,
) -> None:
    """Merge the denoised mates, build the sequence table and remove chimeras.

    R: mergePairs, makeSequenceTable and removeBimeraDenovo. Tables have one row per ASV and
    one column per sample (sorted by name); ASV ids are assigned in --out-seqtab and kept in
    --out-seqtab-nochim. Prints one TSV row per sample: sample, merged, nonchim. Exits with
    code 3 when no sample has any merged read (the tables are still written).

    \b
    Example:
      pydada2 merge --in-fwd asv/*_F.dada.json.gz --in-rev asv/*_R.dada.json.gz \\
        --out-seqtab seqtab.tsv --out-seqtab-nochim seqtab_nochim.tsv --out-fasta asvs.fa
    """  # noqa: D301 - click: \b keeps line breaks; \\ shows a backslash
    with _session(verbose):
        _prepare_outputs(
            [*in_fwd, *in_rev], [out_seqtab, out_seqtab_nochim, out_fasta, out_mergers, stats]
        )
        fwd = _load_mates(in_fwd, "fwd", "--in-fwd")
        rev = _load_mates(in_rev, "rev", "--in-rev")
        if set(fwd) != set(rev):
            missing = sorted(set(fwd) ^ set(rev))
            raise click.UsageError(
                f"Samples without both mates: {', '.join(missing)} (pairing is by the sample "
                "name stored in the .dada.json files)."
            )
        names = sorted(fwd)
        mergers = merge_pairs(
            [fwd[n].result for n in names],
            [fwd[n].derep for n in names],
            [rev[n].result for n in names],
            [rev[n].derep for n in names],
            min_overlap=min_overlap,
            max_mismatch=max_mismatch,
            return_rejects=return_rejects,
            just_concatenate=just_concatenate,
            trim_overhang=trim_overhang,
            n_jobs=threads,
            verbose=verbose,
        )
        assert isinstance(mergers, list)
        accepted = {n: m.loc[m["accept"].astype(bool)] for n, m in zip(names, mergers, strict=True)}
        seqtab = make_sequence_table(
            accepted,
            order_by=None if order_by == "none" else order_by,  # type: ignore[arg-type]
        )
        ids = {str(s): f"ASV{i + 1}" for i, s in enumerate(seqtab.columns)}
        if seqtab.shape[1] > 0:
            nochim = BimeraRemover(
                method=chimera_method,  # type: ignore[arg-type]
                min_fold_parent_over_abundance=min_fold_parent_over_abundance,
                min_parent_abundance=min_parent_abundance,
                min_sample_fraction=min_sample_fraction,
                ignore_n_negatives=ignore_n_negatives,
                allow_one_off=allow_one_off,
                min_one_off_parent_distance=min_one_off_parent_distance,
                max_shift=max_shift,
                n_jobs=threads,
                verbose=verbose,
            ).fit_transform(seqtab)
        else:
            nochim = seqtab
        write_seqtab(seqtab, out_seqtab, asv_ids=list(ids.values()))
        write_seqtab(nochim, out_seqtab_nochim, asv_ids=[ids[str(s)] for s in nochim.columns])
        if out_fasta is not None:
            write_fasta([(ids[str(s)], str(s)) for s in nochim.columns], out_fasta)
        if out_mergers is not None:
            frames = []
            for n, m in zip(names, mergers, strict=True):
                m = m.copy()
                m.insert(0, "asv_id", [ids.get(str(s), "") for s in m["sequence"]])
                m.insert(0, "sample", n)
                frames.append(m)
            _write_tsv(pd.concat(frames, ignore_index=True), out_mergers)
        if verbose:
            logger.info(
                "%d of %d ASVs (%.1f%% of the reads) are not chimeric.",
                nochim.shape[1],
                seqtab.shape[1],
                100 * nochim.to_numpy().sum() / max(1, seqtab.to_numpy().sum()),
            )
        rows = [
            {
                "sample": n,
                "merged": int(seqtab.loc[n].to_numpy().sum()),
                "nonchim": int(nochim.loc[n].to_numpy().sum()),
            }
            for n in names
        ]
        _emit_stats(rows, ["sample", "merged", "nonchim"], stats)
        if not any(r["merged"] for r in rows):
            logger.warning("No sample has any merged read.")
            _no_data()


# --------------------------------------------------------------------- assign-taxonomy ---
def _load_model(path: Path) -> TaxonomyClassifier:
    with open(path, "rb") as fh:
        model = pickle.load(fh)  # trusted: written by --save-model
    if not isinstance(model, TaxonomyClassifier) or not hasattr(model, "model_"):
        raise click.BadParameter(
            f"{path} is not a fitted TaxonomyClassifier (written by --save-model).",
            param_hint="--ref",
        )
    return model


@cli.command("assign-taxonomy")
@click.option(
    "--in",
    "in_",
    type=_IN_FILE,
    required=True,
    help="Sequence table TSV (from merge; only its asv_id and sequence columns are used).",
)
@click.option(
    "--ref",
    type=_IN_FILE,
    required=True,
    help="DADA2-formatted training FASTA (ids are ';'-separated taxonomies, e.g. "
    "silva_nr99_v138.2_toGenus_trainset.fa.gz), or a fitted model saved with "
    "--save-model (a .pkl file).",
)
@click.option(
    "--save-model",
    type=_OUT_FILE,
    default=None,
    help="Pickle the fitted classifier, to skip the training next time (~1.3 GB for SILVA).",
)
@click.option(
    "--out-taxa",
    type=_OUT_FILE,
    required=True,
    help="Taxonomy table (TSV: asv_id, sequence, one column per rank; empty when unassigned).",
)
@click.option(
    "--out-bootstrap",
    type=_OUT_FILE,
    default=None,
    help="Bootstrap confidences (0-100) of the best hit at every rank (same layout).",
)
@click.option(
    "--min-boot",
    type=click.IntRange(0, 100),
    default=50,
    show_default=True,
    help="Minimum bootstrap confidence to assign a rank.",
)
@click.option(
    "--try-rc/--no-try-rc",
    default=False,
    show_default=True,
    help="Also try the reverse complement of each sequence and keep the better "
    "orientation (also used for --species-ref).",
)
@click.option(
    "--tax-levels",
    cls=GreedyOption,
    metavar="L...",
    default=DEFAULT_TAX_LEVELS,
    help="Rank names (column names of the output), from the top level down.  [default: "
    + " ".join(DEFAULT_TAX_LEVELS)
    + "]",
)
@click.option(
    "--seed",
    type=int,
    default=0,
    show_default=True,
    help="Random seed of the bootstrap. Differs from R (which has no seed): a fixed "
    "default keeps reruns reproducible.",
)
@click.option(
    "--species-ref",
    type=_IN_FILE,
    default=None,
    help="Species reference FASTA (ids '<id> <Genus> <species>', e.g. "
    "silva_v138.2_assignSpecies.fa.gz): adds a Species column by exact matching where "
    "the genus agrees (R: addSpecies).",
)
@click.option(
    "--allow-multiple",
    type=_AllowMultiple(),
    is_flag=False,
    flag_value="all",
    default=None,
    metavar="[N]",
    help="Report up to N equally matching species, separated by '/' (--allow-multiple "
    "alone: all). Default: only unambiguous species.",
)
@_stats_option
@_common_options
def assign_taxonomy_cmd(
    in_: Path,
    ref: Path,
    save_model: Path | None,
    out_taxa: Path,
    out_bootstrap: Path | None,
    min_boot: int,
    try_rc: bool,
    tax_levels: tuple[str, ...],
    seed: int,
    species_ref: Path | None,
    allow_multiple: bool | int | None,
    stats: Path | None,
    threads: int,
    verbose: int,
) -> None:
    """Assign taxonomy to the ASVs of a sequence table (R: assignTaxonomy, addSpecies).

    Naive Bayesian classifier (Wang et al. 2007) with bootstrap confidences, plus exact
    species matching with --species-ref. Rows keep the order of the input table. Prints one
    TSV row per rank: rank, assigned, total. Exits with code 3 if the table has no ASVs.

    \b
    Example:
      pydada2 assign-taxonomy --in seqtab_nochim.tsv --ref silva_train.fa.gz \\
        --out-taxa taxa.tsv --out-bootstrap bootstrap.tsv --species-ref silva_species.fa.gz
    """  # noqa: D301 - click: \b keeps line breaks; \\ shows a backslash
    with _session(verbose):
        _prepare_outputs(
            [in_, ref, *([species_ref] if species_ref else [])],
            [save_model, out_taxa, out_bootstrap, stats],
        )
        table = pd.read_csv(in_, sep="\t", dtype=str, keep_default_na=False)
        if not {"asv_id", "sequence"} <= set(table.columns):
            raise click.BadParameter(
                "the table needs 'asv_id' and 'sequence' columns.", param_hint="--in"
            )
        ids, seqs = table["asv_id"].tolist(), table["sequence"].tolist()
        settings: dict[str, Any] = {
            "min_boot": min_boot,
            "try_rc": try_rc,
            "tax_levels": tuple(tax_levels),
            "random_state": seed,
            "n_jobs": threads,
            "verbose": verbose,
        }
        if ref.name.endswith(".pkl"):
            clf = _load_model(ref).set_params(**settings)
        else:
            clf = TaxonomyClassifier(**settings).fit(ref)
        if save_model is not None:
            with open(save_model, "wb") as fh:
                pickle.dump(clf, fh, protocol=pickle.HIGHEST_PROTOCOL)
        if seqs:
            taxa, boot = clf.predict(seqs, return_bootstrap=True)
            if species_ref is not None:
                taxa = add_species(
                    taxa,
                    species_ref,
                    allow_multiple=False if allow_multiple is None else allow_multiple,
                    try_rc=try_rc,
                    n_jobs=threads,
                    verbose=verbose,
                )
        else:
            cols = clf._columns() + (["Species"] if species_ref is not None else [])
            taxa = pd.DataFrame(columns=list(dict.fromkeys(cols)), dtype=object)
            boot = pd.DataFrame(columns=clf._columns(), dtype=np.int64)
        for frame, path in ((taxa, out_taxa), (boot, out_bootstrap)):
            if path is None:
                continue
            out = frame.reset_index(drop=True)
            out.insert(0, "sequence", seqs)
            out.insert(0, "asv_id", ids)
            _write_tsv(out, path)
        rows = [
            {"rank": str(c), "assigned": int(taxa[c].notna().sum()), "total": len(seqs)}
            for c in taxa.columns
        ]
        _emit_stats(rows, ["rank", "assigned", "total"], stats)
        if not seqs:
            logger.warning("The sequence table has no ASVs.")
            _no_data()


def main(args: Sequence[str] | None = None) -> None:
    """Run the ``pydada2`` command line (the console-script entry point).

    Parameters
    ----------
    args : sequence of str, optional
        Command-line arguments (default: ``sys.argv[1:]``).
    """
    cli.main(args=list(args) if args is not None else None, prog_name="pydada2")
