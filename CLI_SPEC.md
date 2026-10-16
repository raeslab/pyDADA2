# CLI

A `click`-based CLI that runs the DADA2 paired-end workflow as separate steps, orchestrated (and
parallelised) by Nextflow. Decided in interviews on 2026-10-16.

```bash
pydada2 trim-and-filter  --in-fwd F.fq.gz --in-rev R.fq.gz --out-fwd filt_F.fq.gz --out-rev filt_R.fq.gz [opts]  # per sample
pydada2 learn-errors     --in-fwd filt/*_F.fq.gz --in-rev filt/*_R.fq.gz --out errors.pkl [opts]                  # per batch
pydada2 infer-asv        --in-fwd filt_F.fq.gz --in-rev filt_R.fq.gz --errors errors.pkl \
                         --out-fwd sam1_F.dada.json.gz --out-rev sam1_R.dada.json.gz [opts]                       # per sample
pydada2 merge            --in-fwd asv/*_F.dada.json.gz --in-rev asv/*_R.dada.json.gz \
                         --out-seqtab seqtab.tsv --out-seqtab-nochim seqtab_nochim.tsv [opts]                     # per batch, also removes chimeras
pydada2 assign-taxonomy  --in seqtab_nochim.tsv --ref silva_train.fa.gz --out-taxa taxa.tsv [opts]                # once
```

## Conventions (all subcommands)

- **git-style:** a `pydada2` console script with subcommands. `python -m pydada2 ...` is
  equivalent. Subcommands and options are **kebab-case**, and a Python argument maps to its
  kebab-case name (`trunc_len` → `--trunc-len`, `min_fold_parent_over_abundance` →
  `--min-fold-parent-over-abundance`).
- **Named options only.** Every input and output is an option. There are no positional
  arguments, and order is free.
- **Shared names:** `--in-fwd`/`--in-rev` and `--out-fwd`/`--out-rev` for per-mate files,
  `--in`/`--out` for single files, `--out-<what>` for further outputs.
- **File lists** (batch steps) take one or more space-separated paths after a single flag, so
  shell globs work (`--in-fwd filt/*_F.fq.gz`). The list ends at the next `-`-prefixed token.
  click has no variadic options, so this needs a small greedy `Option` subclass.
- **Per-mate parameters** take `FWD [REV]`: one value applies to both mates, two set them
  separately (R recycling). `inf`/`-inf` are accepted where the Python API accepts them.
- **Defaults** are R's (as in the Python API) unless noted. Every Python parameter is exposed.
- **Sample names:** `--sample-name NAME` overrides, else the name is derived from the
  `--in-fwd` basename: strip `.gz`, then `.fastq`/`.fq`/`.dada.json`, then a trailing
  `_R1`/`_1`/`_F` (`sam1_R1.fastq.gz` → `sam1`). Batch steps take the name stored in each
  input file.
- **Output files** are gzipped when the path ends in `.gz` (no `--compress` flag). Missing
  directories are created, and outputs must differ from inputs.
- **stdout carries only machine-readable TSV** (read tracking, header line + rows). `--stats
  PATH` also writes the same TSV to a file (Nextflow: `collectFile(keepHeader: true)`).
  Logs and progress go to stderr (`-v/--verbose`, repeatable).
- **Threads:** `--threads N`, default **1** (Nextflow passes `${task.cpus}`).
- **Random seeds:** `--seed N` where the step is stochastic.
- **Exit codes:** `0` ok, `1` runtime error, `2` usage error (click), `3` no data (stats are still
  written). In Nextflow: `errorStrategy { task.exitStatus == 3 ? 'ignore' : 'terminate' }`.
- **Tables** are TSV. ASV tables are keyed by `asv_id` + `sequence` (rows = ASVs).

## trim-and-filter (per sample)

Wraps `filter_and_trim(in_fwd, out_fwd, in_rev, out_rev, ...)` for one read pair.

```bash
pydada2 trim-and-filter \
  --in-fwd sam1_R1.fastq.gz --in-rev sam1_R2.fastq.gz \
  --out-fwd filtered/sam1_F.fq.gz --out-rev filtered/sam1_R.fq.gz \
  --trunc-len 240 160 --max-ee 2 --trunc-q 2 \
  --sample-name sam1 --stats sam1.filter_stats.tsv --threads 2
```

| Option | Default | Description |
|---|---|---|
| `--in-fwd PATH` | required | Forward-read FASTQ (optionally gzipped). |
| `--in-rev PATH` | required | Matching reverse-read FASTQ. |
| `--out-fwd PATH` | required | Filtered forward output. |
| `--out-rev PATH` | required | Filtered reverse output. |

Per-mate (`FWD [REV]`):

| Option | Default | Description |
|---|---|---|
| `--trunc-q` | `2` | Truncate at the first base with quality <= this. |
| `--trunc-len` | `0` | Truncate after N bases; shorter reads are discarded (0 = off). |
| `--trim-left` | `0` | Bases removed from the start. |
| `--trim-right` | `0` | Bases removed from the end. |
| `--max-len` | `inf` | Discard reads longer than this (before trimming). |
| `--min-len` | `20` | Discard reads shorter than this (after trimming). |
| `--max-n` | `0` | Max number of non-ACGT bases. |
| `--min-q` | `0` | Discard reads with a quality <= this (only when > `--trunc-q`, as in R). |
| `--max-ee` | `inf` | Max expected errors. |
| `--rm-lowcomplex` | `0` | Discard reads with sequence complexity below this (0 = off). |

Other:

| Option | Default | Description |
|---|---|---|
| `--rm-phix / --keep-phix` | rm | Remove reads matching phiX. |
| `--orient-fwd SEQ` | none | Keep only pairs starting with SEQ; swap mates when needed. |
| `--match-ids / --no-match-ids` | no | Keep only pairs whose ids match. |
| `--id-sep REGEX` | `\s` | Id separator for `--match-ids`. |
| `--id-field N` | auto | Id field for `--match-ids` (auto-detects CASAVA). |
| `--quality-type` | `auto` | `auto`, `FastqQuality` or `SFastqQuality`. |
| `--chunk-size N` | `100000` | Reads per chunk (Python `n`). |
| `--sample-name NAME` | derived | Name in the stats row. |
| `--stats PATH` | none | Also write the stats TSV here. |
| `--threads N` | `1` | Worker threads. |
| `-v, --verbose` | off | Log to stderr. |

stdout:

```
sample	reads_in	reads_out
sam1	7113	6996
```

**Zero reads pass:** the stats row is still written (`reads_out` = 0), no FASTQ is written, and
the command exits with code **3**.

## learn-errors (per batch)

Wraps `learn_errors(fwd_files)` and `learn_errors(rev_files)`.

```bash
pydada2 learn-errors \
  --in-fwd filtered/*_F.fq.gz --in-rev filtered/*_R.fq.gz \
  --out errors.pkl --out-tsv errors.tsv.gz --plot errors.pdf --threads 8
```

| Option | Default | Description |
|---|---|---|
| `--in-fwd PATH...` | required | Filtered forward FASTQs. |
| `--in-rev PATH...` | required | Filtered reverse FASTQs (same count; order doesn't matter). |
| `--out PATH` | required | Pickle of `{"err_f": ndarray, "err_r": ndarray}` (the 16 × Q `err_` matrices). |
| `--out-tsv PATH` | none | Long-format dump to read: `mate transition qual rate` (`A2C` ...). |
| `--plot PATH` | none | `plot_errors` figure for both mates (format from extension: png/pdf/svg). |
| `--n-bases N` | `1e8` | Min number of bases to learn from (per mate). |
| `--error-model` | `loess` | `loess`, `noqual`, `binned` or `pacbio`. |
| `--quality-bins Q...` | none | Bins for `--error-model binned`. |
| `--randomize / --no-randomize` | no | Pick samples in random order. |
| `--seed N` | none | Seed for `--randomize`. |
| `--max-consist N` | `10` | Max self-consistency iterations. |
| `--omega-c X` | `0` | |
| DADA options | | All `DadaParams` options, as in `infer-asv` (`--band-size`, ...). |
| `--quality-type` | `auto` | |
| `--threads N` | `1` | |
| `-v, --verbose` | off | |

stdout (one row per mate):

```
mate	n_samples	n_reads	n_bases	n_iter	converged
fwd	2	6996	1679040	4	true
rev	2	6996	1119360	5	true
```

A non-converged run is a stderr warning (exit 0), not an error.

## infer-asv (per sample)

Dereplicates both mates and runs `DADA(err=...)` on each (the error pickle from `learn-errors`).
It does **not** merge. The outputs keep everything `merge_pairs` needs.

```bash
pydada2 infer-asv \
  --in-fwd filtered/sam1_F.fq.gz --in-rev filtered/sam1_R.fq.gz --errors errors.pkl \
  --out-fwd asv/sam1_F.dada.json.gz --out-rev asv/sam1_R.dada.json.gz \
  --sample-name sam1 --threads 2
```

| Option | Default | Description |
|---|---|---|
| `--in-fwd PATH` | required | Filtered forward FASTQ. |
| `--in-rev PATH` | required | Filtered reverse FASTQ. |
| `--errors PATH` | required | Error pickle from `learn-errors`. |
| `--out-fwd PATH` | required | Forward result (`.dada.json.gz`). |
| `--out-rev PATH` | required | Reverse result (`.dada.json.gz`). |
| `--sample-name NAME` | derived | Stored in the outputs. |
| `--priors-fwd/--priors-rev FASTA` | none | Prior sequences. |
| `--chunk-size N` | `1000000` | derep chunk size (Python `n`). |
| `--quality-type` | `auto` | |
| DADA options | | One per `DadaParams` field, R defaults: `--omega-a 1e-40`, `--omega-p 1e-4`, `--omega-c 1e-40`, `--detect-singletons/--no-...`, `--use-kmers/--no-...`, `--kdist-cutoff 0.42`, `--match 5`, `--mismatch -4`, `--gap-penalty -8`, `--homopolymer-gap-penalty`, `--band-size 16`, `--max-clust 0`, `--min-fold 1`, `--min-hamming 1`, `--min-abundance 1`, `--use-quals/--no-...`, `--gapless/--no-...`, `--greedy/--no-...`, `--vectorized-alignment/--no-...`, `--sse 2`. |
| `--stats PATH` | none | |
| `--threads N` | `1` | |
| `-v, --verbose` | off | |

Pooling (`pool=True`/`"pseudo"`) is out of scope: it is not per-sample.

**Output format `pydada2.dada/1`** (gzipped JSON, one file per mate; inspect with
`zcat x.dada.json.gz | jq .clusters`), read back with a new public `pydada2.read_dada_json()`
(and written by `write_dada_json()`):

```json
{
  "format": "pydada2.dada/1",
  "pydada2_version": "x.y.z",
  "sample": "sam1",
  "mate": "fwd",
  "source": "sam1_F.fq.gz",
  "clusters": {"sequence": [...], "abundance": [...], "n0": [...], "pval": [...], "birth_from": [...], ...},
  "dada_map": [0, 0, 1, -1, ...],
  "derep": {"sequence": [...], "abundance": [...], "map": [...]},
  "params": {"omega_a": 1e-40, ...}
}
```

Columns are stored as column lists (compact), `NaN`/`inf` as `null`. Large per-position arrays
(`quality`, derep `quals`) aren't needed for merging and are left out (`--full` could add them
later).

stdout:

```
sample	reads	uniques_fwd	uniques_rev	denoised_fwd	denoised_rev
sam1	6996	1979	1660	6865	6855
```

Zero input reads: both outputs are written with empty tables, stats are printed, exit code **3**.

## merge (per batch; also removes chimeras)

Loads the `.dada.json.gz` files, runs `merge_pairs` per sample, then `make_sequence_table` and
`BimeraRemover`.

```bash
pydada2 merge \
  --in-fwd asv/*_F.dada.json.gz --in-rev asv/*_R.dada.json.gz \
  --out-seqtab seqtab.tsv --out-seqtab-nochim seqtab_nochim.tsv \
  --out-fasta asvs_nochim.fa --min-overlap 12 --threads 8
```

| Option | Default | Description |
|---|---|---|
| `--in-fwd PATH...` | required | Forward `.dada.json.gz` files. |
| `--in-rev PATH...` | required | Reverse files. Paired by the stored `sample` (not position); error on a missing/duplicate mate. |
| `--out-seqtab PATH` | required | Sequence table before chimera removal. |
| `--out-seqtab-nochim PATH` | required | Sequence table after chimera removal. |
| `--out-fasta PATH` | none | Non-chimeric ASVs as FASTA (`>ASV1`). |
| `--out-mergers PATH` | none | Long table of every merger per sample (incl. rejects when `--return-rejects`). |
| `--min-overlap N` | `12` | |
| `--max-mismatch N` | `0` | |
| `--return-rejects / --no-...` | no | Keep rejected pairs in `--out-mergers`. |
| `--just-concatenate / --no-...` | no | Concatenate with 10 Ns instead of merging. |
| `--trim-overhang / --no-...` | no | |
| `--order-by` | `abundance` | `abundance`, `nsamples` or `none` (row order and ASV numbering). |
| `--chimera-method` | `consensus` | `consensus`, `pooled` or `per-sample`. |
| `--min-fold-parent-over-abundance X` | 1.5 / 2 | Default per method (as in R). |
| `--min-parent-abundance X` | 2 / 8 | Default per method. |
| `--min-sample-fraction X` | `0.9` | consensus only. |
| `--ignore-n-negatives N` | `1` | |
| `--allow-one-off / --no-...` | no | |
| `--min-one-off-parent-distance N` | `4` | |
| `--max-shift N` | `16` | |
| `--stats PATH` | none | |
| `--threads N` | `1` | |
| `-v, --verbose` | off | |

**Seqtab format** (both tables): TSV, one row per ASV, one column per sample (sorted by name):

```
asv_id	sequence	sam1	sam2
ASV1	TACGGAGG...	579	345
ASV2	TACGGAGG...	470	352
```

- This is the transpose of R's `seqtab`, chosen because it is readable and `cut`/`awk`-friendly
  (no 250 bp column headers). Tabs can't occur in sequences.
- `asv_id`s are assigned once in `seqtab` and **kept** in `seqtab_nochim` (removed chimeras
  leave gaps), so the two tables join.
- The reader `pydada2.read_seqtab()` returns R's orientation (samples × sequences).

stdout:

```
sample	merged	nonchim
sam1	6540	6528
```

Exit code **3** when no sample has any merged read.

## assign-taxonomy (once, on the non-chimeric table)

Wraps `TaxonomyClassifier(...).fit(ref).predict(seqs, return_bootstrap=True)`, plus
`add_species` when `--species-ref` is given.

```bash
pydada2 assign-taxonomy \
  --in seqtab_nochim.tsv --ref silva_train.fa.gz \
  --out-taxa taxa.tsv --out-bootstrap bootstrap.tsv \
  --species-ref silva_species.fa.gz --min-boot 50 --seed 0 --threads 8
```

| Option | Default | Description |
|---|---|---|
| `--in PATH` | required | Seqtab TSV (uses the `asv_id` and `sequence` columns). |
| `--ref PATH` | required | DADA2-formatted training FASTA, **or** a fitted model `.pkl` (detected by extension). |
| `--save-model PATH` | none | Pickle the fitted model for reuse (~1.3 GB for SILVA). |
| `--out-taxa PATH` | required | Taxonomy table. |
| `--out-bootstrap PATH` | none | Bootstrap values (same shape). |
| `--min-boot N` | `50` | Min bootstrap to assign a rank. |
| `--try-rc / --no-try-rc` | no | Also try the reverse complement. |
| `--tax-levels L...` | `Kingdom Phylum Class Order Family Genus Species` | Rank names. |
| `--seed N` | `0` | **Differs from R** (no seed): a fixed default keeps Nextflow runs reproducible. |
| `--species-ref PATH` | none | Exact-match species FASTA → adds/fills `Species` (`add_species`). |
| `--allow-multiple N` | off | Report up to N species (`--allow-multiple` alone: all). |
| `--threads N` | `1` | |
| `-v, --verbose` | off | |

Table format (rows in the input's order; an unassigned rank is an empty cell):

```
asv_id	sequence	Kingdom	Phylum	Class	Order	Family	Genus	Species
ASV1	TACG...	Bacteria	Bacteroidota	Bacteroidia	Bacteroidales	Muribaculaceae
ASV2	TACG...	Bacteria	Firmicutes	Clostridia	Lachnospirales	Lachnospiraceae	Lachnoclostridium
```

`bootstrap.tsv` has the same columns with integers 0–100 (no Species column; exact matching has
no bootstrap). stdout is a per-rank summary:

```
rank	assigned	total
Kingdom	240	240
Genus	180	240
```

## Still to do

- A CLI page in the documentation: how to run each step, a pipeline diagram, a table of the
  files each step expects and produces, and a minimal Nextflow example.
- Pending decisions: none.
