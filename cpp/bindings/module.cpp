// nanobind bindings for the pydada2 C++ core (private module pydada2._core).
//
// One function per upstream DADA2 Rcpp export (snake_case), plus TaxonomyModel and a few
// private test helpers. The public, documented API lives in Python (src/pydada2).
// Heavy calls release the GIL; log/interrupt callbacks re-acquire it (main thread only).

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <atomic>
#include <climits>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "api.h"
#include "parallel.h"
#include "ported/ported.h"
#include "runtime.h"
#include "stats.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

constexpr const char *CORE_VERSION = "0.1.0";
constexpr const char *UPSTREAM_VERSION = "1.36.0";

template <typename T>
using CArray = nb::ndarray<const T, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

// Integer matrices for the R-integer kernels (table_bimera, matrix_ee). Each is bound twice: an
// int32 C-contiguous array is used in place; any other numeric array is converted to int64 by
// nanobind (that overload is registered first, so implicit conversions go there) and copied
// with a range check, instead of nanobind's narrowing cast silently wrapping (2**32 + 5 -> 5).
using I64Array = nb::ndarray<const int64_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

std::vector<int32_t> checked_int32(const I64Array &a, const char *name) {
  std::vector<int32_t> out(a.size());
  const int64_t *p = a.data();
  for (size_t i = 0; i < out.size(); ++i) {
    if (p[i] < INT32_MIN || p[i] > INT32_MAX)
      throw std::invalid_argument(std::string(name) + " has values outside the int32 range.");
    out[i] = static_cast<int32_t>(p[i]);
  }
  return out;
}

// Move a std::vector into a NumPy array (1-D or 2-D) without copying.
template <typename T>
nb::ndarray<nb::numpy, T> to_numpy(std::vector<T> &&v, std::initializer_list<size_t> shape) {
  auto *heap = new std::vector<T>(std::move(v));
  nb::capsule owner(heap, [](void *p) noexcept { delete static_cast<std::vector<T> *>(p); });
  return nb::ndarray<nb::numpy, T>(heap->data(), shape.size(), shape.begin(), owner);
}

template <typename T>
nb::ndarray<nb::numpy, T> to_numpy(std::vector<T> &&v) {
  size_t n = v.size();
  return to_numpy(std::move(v), {n});
}

nb::ndarray<nb::numpy, bool> bools_to_numpy(const std::vector<bool> &v) {
  std::vector<bool> copy = v;
  auto *heap = new std::vector<uint8_t>(copy.begin(), copy.end());
  nb::capsule owner(heap, [](void *p) noexcept { delete static_cast<std::vector<uint8_t> *>(p); });
  size_t shape[1] = {heap->size()};
  return nb::ndarray<nb::numpy, bool>(reinterpret_cast<bool *>(heap->data()), 1, shape, owner);
}

// Installs log/interrupt handlers that call back into Python with the GIL held.
struct PyHandlers {
  dada2::ScopedHandlers scoped;
  explicit PyHandlers(std::optional<nb::callable> log)
      : scoped(make_log(std::move(log)), make_interrupt()) {}

  static dada2::LogHandler make_log(std::optional<nb::callable> log) {
    if (!log) return nullptr;
    auto fn = std::make_shared<nb::callable>(std::move(*log));
    return [fn](const std::string &msg) {
      nb::gil_scoped_acquire gil;
      (*fn)(msg);
    };
  }
  static dada2::InterruptHandler make_interrupt() {
    return []() {
      nb::gil_scoped_acquire gil;
      return PyErr_CheckSignals() != 0;  // leaves KeyboardInterrupt (or other) set
    };
  }
};

// Dereplicator with a per-instance mutex: add_chunk/finish release the GIL, so concurrent calls
// on one instance from different Python threads must not race on its hash table.
struct LockedDereplicator {
  dada2::Dereplicator impl;
  std::mutex mutex;
};

nb::ndarray<nb::numpy, double> matrix_ee_py(const int32_t *quals, size_t nrow, size_t ncol) {
  std::vector<double> out;
  {
    nb::gil_scoped_release release;
    out = dada2::matrix_ee(quals, static_cast<int>(nrow), static_cast<int>(ncol));
  }
  return to_numpy(std::move(out));
}

std::pair<nb::ndarray<nb::numpy, int>, nb::ndarray<nb::numpy, int>> table_bimera_py(
    const int32_t *mat, size_t nrow, size_t ncol, const std::vector<std::string> &seqs,
    double min_fold, int min_abund, bool allow_one_off, int min_one_off_par_dist, int match,
    int mismatch, int gap_p, int max_shift, int n_threads) {
  std::pair<std::vector<int>, std::vector<int>> r;
  {
    nb::gil_scoped_release release;
    r = dada2::table_bimera(mat, static_cast<int>(nrow), static_cast<int>(ncol), seqs, min_fold,
                            min_abund, allow_one_off, min_one_off_par_dist, match, mismatch, gap_p,
                            max_shift, n_threads);
  }
  return std::make_pair(to_numpy(std::move(r.first)), to_numpy(std::move(r.second)));
}

nb::dict dada_output_to_dict(dada2::DadaOutput &&o) {
  nb::dict d;
  const size_t nclust = o.sequence.size();
  d["sequence"] = nb::cast(std::move(o.sequence));
  d["abundance"] = to_numpy(std::move(o.abundance));
  d["n0"] = to_numpy(std::move(o.n0));
  d["n1"] = to_numpy(std::move(o.n1));
  d["nunq"] = to_numpy(std::move(o.nunq));
  d["pval"] = to_numpy(std::move(o.pval));
  d["birth_from"] = to_numpy(std::move(o.birth_from));
  d["birth_pval"] = to_numpy(std::move(o.birth_pval));
  d["birth_fold"] = to_numpy(std::move(o.birth_fold));
  d["birth_ham"] = to_numpy(std::move(o.birth_ham));
  d["birth_qave"] = to_numpy(std::move(o.birth_qave));
  d["bs_pos"] = to_numpy(std::move(o.bs_pos));
  d["bs_ref"] = nb::cast(std::move(o.bs_ref));
  d["bs_sub"] = nb::cast(std::move(o.bs_sub));
  d["bs_qual"] = to_numpy(std::move(o.bs_qual));
  d["bs_clust"] = to_numpy(std::move(o.bs_clust));
  const size_t tcol = static_cast<size_t>(o.trans_ncol);
  d["trans"] = to_numpy(std::move(o.trans), {16, tcol});
  const size_t qcol = static_cast<size_t>(o.quals_ncol);
  d["cluster_quals"] = to_numpy(std::move(o.cluster_quals), {nclust, qcol});
  d["map"] = to_numpy(std::move(o.map));
  d["pval_raw"] = to_numpy(std::move(o.pval_raw));
  return d;
}

}  // namespace

NB_MODULE(_core, m) {
  m.doc() = "Private C++ core of pydada2 (ported from DADA2). Not a public API.";

  nb::register_exception_translator([](const std::exception_ptr &p, void *) {
    try {
      std::rethrow_exception(p);
    } catch (const dada2::Interrupted &) {
      if (!PyErr_Occurred()) PyErr_SetNone(PyExc_KeyboardInterrupt);
    }
  });

  m.def("core_version", []() { return std::string(CORE_VERSION); },
        "Version of the pydada2 C++ core.");
  m.def("upstream_version", []() { return std::string(UPSTREAM_VERSION); },
        "Version of upstream DADA2 the core was ported from.");

  // ------------------------------------------------------------- stats / runtime helpers ---
  m.def("ppois_upper", &ppois_upper, "n"_a, "lam"_a,
        "Upper Poisson tail P(X > n), X ~ Poisson(lam) (R: ppois(n, lam, lower.tail=FALSE)).");

  m.def(
      "_parallel_selftest",
      [](size_t n, size_t grain, int n_threads, long long fail_at) {
        std::vector<int64_t> out(n, 0);
        std::atomic<int> calls{0};
        {
          nb::gil_scoped_release release;
          dada2::parallel_for(
              0, n,
              [&](size_t b, size_t e) {
                calls++;
                for (size_t i = b; i < e; ++i) {
                  if (fail_at >= 0 && static_cast<long long>(i) == fail_at)
                    throw dada2::Error("parallel failure at " + std::to_string(i));
                  out[i] = static_cast<int64_t>(i) * static_cast<int64_t>(i);
                }
              },
              grain, n_threads);
        }
        return std::make_pair(to_numpy(std::move(out)), calls.load());
      },
      "n"_a, "grain"_a, "n_threads"_a, "fail_at"_a = -1,
      "Test helper: squares 0..n-1 with parallel_for; returns (values, number of chunks).");

  m.def(
      "_interrupt_selftest",
      [](long long n_checks) {
        PyHandlers handlers(std::nullopt);
        nb::gil_scoped_release release;
        for (long long i = 0; i < n_checks; ++i) dada2_check_interrupt();
      },
      "n_checks"_a, "Test helper: calls the interrupt check n_checks times (GIL released).");

  // ------------------------------------------------------------------------ dada_uniques ---
  m.def(
      "dada_uniques",
      [](const std::vector<std::string> &seqs, const std::vector<int> &abundances,
         const std::vector<bool> &priors, CArray<double> err,
         std::optional<CArray<double>> quals, int match, int mismatch, int gap, bool use_kmers,
         double kdist_cutoff, int band_size, double omega_a, double omega_p, double omega_c,
         bool detect_singletons, int max_clust, double min_fold, int min_hamming, int min_abund,
         bool use_quals, bool final_consensus, bool vectorized_alignment, int homo_gap, int sse,
         bool gapless, bool greedy, int n_threads, bool verbose,
         std::optional<nb::callable> log) {
        if (err.shape(0) != 16) throw std::invalid_argument("Error matrix must have 16 rows.");
        const double *qptr = nullptr;
        int qcol = 0;
        if (quals && quals->shape(0) > 0) {
          if (quals->shape(0) != seqs.size())
            throw std::invalid_argument("quals must have one row per unique sequence.");
          qptr = quals->data();
          qcol = static_cast<int>(quals->shape(1));
        }
        dada2::DadaOptions opt;
        opt.match = match; opt.mismatch = mismatch; opt.gap = gap; opt.use_kmers = use_kmers;
        opt.kdist_cutoff = kdist_cutoff; opt.band_size = band_size; opt.omega_a = omega_a;
        opt.omega_p = omega_p; opt.omega_c = omega_c; opt.detect_singletons = detect_singletons;
        opt.max_clust = max_clust; opt.min_fold = min_fold; opt.min_hamming = min_hamming;
        opt.min_abund = min_abund; opt.use_quals = use_quals; opt.final_consensus = final_consensus;
        opt.vectorized_alignment = vectorized_alignment; opt.homo_gap = homo_gap; opt.sse = sse;
        opt.gapless = gapless; opt.greedy = greedy;
        dada2::DadaOutput out;
        {
          PyHandlers handlers(std::move(log));
          nb::gil_scoped_release release;
          out = dada2::dada_uniques(seqs, abundances, priors, err.data(),
                                    static_cast<int>(err.shape(1)), qptr, qcol, opt, n_threads,
                                    verbose);
        }
        return dada_output_to_dict(std::move(out));
      },
      "seqs"_a, "abundances"_a, "priors"_a, "err"_a, "quals"_a.none(), "match"_a = 5,
      "mismatch"_a = -4, "gap"_a = -8, "use_kmers"_a = true, "kdist_cutoff"_a = 0.42,
      "band_size"_a = 16, "omega_a"_a = 1e-40, "omega_p"_a = 1e-4, "omega_c"_a = 1e-40,
      "detect_singletons"_a = false, "max_clust"_a = 0, "min_fold"_a = 1.0, "min_hamming"_a = 1,
      "min_abund"_a = 1, "use_quals"_a = true, "final_consensus"_a = false,
      "vectorized_alignment"_a = true, "homo_gap"_a = -8, "sse"_a = 2, "gapless"_a = true,
      "greedy"_a = true, "n_threads"_a = 1, "verbose"_a = false, "log"_a = nb::none(),
      "Run the DADA2 sample-inference algorithm on dereplicated sequences.\n\n"
      "err is 16 x Q (rows A2A..T2T); quals is n_uniques x max_len (or None). Returns a dict of\n"
      "arrays (clustering columns, birth_subs columns, trans, cluster_quals, map, pval_raw);\n"
      "indices are 0-based, -1/NaN mean NA.");

  // ---------------------------------------------------------------------------- evaluate ---
  m.def("nwalign", &dada2::nwalign, "s1"_a, "s2"_a, "match"_a, "mismatch"_a, "gap_p"_a,
        "homo_gap_p"_a, "band"_a, "endsfree"_a, nb::call_guard<nb::gil_scoped_release>(),
        "Needleman-Wunsch alignment (ACGT only). Returns the two aligned strings.");
  m.def("nwvec", &dada2::nwvec, "s1"_a, "s2"_a, "match"_a, "mismatch"_a, "gap_p"_a, "band"_a,
        "endsfree"_a, nb::call_guard<nb::gil_scoped_release>(),
        "Vectorized Needleman-Wunsch alignment of s1[i] with s2[i].");
  m.def("eval_pair", &dada2::eval_pair, "s1"_a, "s2"_a,
        "Count (match, mismatch, indel) in the internal part of an alignment.");
  m.def("pair_consensus", &dada2::pair_consensus, "s1"_a, "s2"_a, "prefer"_a, "trim_overhang"_a,
        "Consensus of two aligned sequences (prefer=1|2 wins mismatches), gaps removed.");
  m.def(
      "is_acgt", [](const std::vector<std::string> &seqs) { return bools_to_numpy(dada2::is_acgt(seqs)); },
      "seqs"_a, "Whether each sequence consists only of A/C/G/T.");
  m.def("kmer_dist", &dada2::kmer_dist, "s1"_a, "s2"_a, "kmer_size"_a,
        nb::call_guard<nb::gil_scoped_release>(), "k-mer distance between s1[i] and s2[i].");
  m.def("kord_dist", &dada2::kord_dist, "s1"_a, "s2"_a, "kmer_size"_a, "sse"_a,
        nb::call_guard<nb::gil_scoped_release>(), "Ordered k-mer distance between s1[i] and s2[i].");
  m.def("kmer_matches", &dada2::kmer_matches, "s1"_a, "s2"_a, "kmer_size"_a,
        nb::call_guard<nb::gil_scoped_release>(), "Positional k-mer matches between s1[i] and s2[i].");
  m.def("kdist_matches", &dada2::kdist_matches, "s1"_a, "s2"_a, "kmer_size"_a,
        nb::call_guard<nb::gil_scoped_release>(), "Shared k-mer counts between s1[i] and s2[i].");

  // ------------------------------------------------------------------------------ filter ---
  m.def("match_ref", &dada2::match_ref, "seqs"_a, "ref"_a, "word_size"_a, "non_overlapping"_a,
        nb::call_guard<nb::gil_scoped_release>(),
        "Number of word_size-mers of each sequence found in the (circular) reference.");
  m.attr("NA_QUAL") = dada2::NA_QUAL;
  m.def(
      "matrix_ee",
      [](I64Array quals) {
        const std::vector<int32_t> q = checked_int32(quals, "quals");
        return matrix_ee_py(q.data(), quals.shape(0), quals.shape(1));
      },
      "quals"_a);
  m.def(
      "matrix_ee",
      [](CArray<int32_t> quals) {
        return matrix_ee_py(quals.data(), quals.shape(0), quals.shape(1));
      },
      "quals"_a,
      "Expected errors per row of an int32 quality matrix (NA_QUAL ends a row). Other integer\n"
      "dtypes are accepted if every value fits in int32.");

  // ----------------------------------------------------------------------------- chimera ---
  m.def("is_bimera", &dada2::is_bimera, "sq"_a, "pars"_a, "allow_one_off"_a,
        "min_one_off_par_dist"_a, "match"_a, "mismatch"_a, "gap_p"_a, "max_shift"_a,
        nb::call_guard<nb::gil_scoped_release>(),
        "Whether sq is a bimera of two of the given parents.");
  m.def(
      "table_bimera",
      [](I64Array mat, const std::vector<std::string> &seqs, double min_fold, int min_abund,
         bool allow_one_off, int min_one_off_par_dist, int match, int mismatch, int gap_p,
         int max_shift, int n_threads) {
        const std::vector<int32_t> v = checked_int32(mat, "mat");
        return table_bimera_py(v.data(), mat.shape(0), mat.shape(1), seqs, min_fold, min_abund,
                               allow_one_off, min_one_off_par_dist, match, mismatch, gap_p,
                               max_shift, n_threads);
      },
      "mat"_a, "seqs"_a, "min_fold"_a, "min_abund"_a, "allow_one_off"_a,
      "min_one_off_par_dist"_a, "match"_a, "mismatch"_a, "gap_p"_a, "max_shift"_a,
      "n_threads"_a = 1);
  m.def(
      "table_bimera",
      [](CArray<int32_t> mat, const std::vector<std::string> &seqs, double min_fold, int min_abund,
         bool allow_one_off, int min_one_off_par_dist, int match, int mismatch, int gap_p,
         int max_shift, int n_threads) {
        return table_bimera_py(mat.data(), mat.shape(0), mat.shape(1), seqs, min_fold, min_abund,
                               allow_one_off, min_one_off_par_dist, match, mismatch, gap_p,
                               max_shift, n_threads);
      },
      "mat"_a, "seqs"_a, "min_fold"_a, "min_abund"_a, "allow_one_off"_a,
      "min_one_off_par_dist"_a, "match"_a, "mismatch"_a, "gap_p"_a, "max_shift"_a,
      "n_threads"_a = 1,
      "Per-sequence (nflag, nsam) bimera votes for a samples x sequences int32 table. Other\n"
      "integer dtypes are accepted if every value fits in int32.");

  // ---------------------------------------------------------------------------- taxonomy ---
  nb::class_<dada2::TaxonomyModel>(m, "TaxonomyModel",
                                   "Reference k-mer model for the naive Bayesian classifier.")
      .def(
          "__init__",
          [](dada2::TaxonomyModel *self, const std::vector<std::string> &refs,
             const std::vector<int> &ref_to_genus, CArray<int32_t> genusmat) {
            std::vector<int> gm(genusmat.data(), genusmat.data() + genusmat.size());
            const int ngenus = static_cast<int>(genusmat.shape(0));
            const int nlevel = static_cast<int>(genusmat.shape(1));
            nb::gil_scoped_release release;
            new (self) dada2::TaxonomyModel(refs, ref_to_genus, gm, ngenus, nlevel);
          },
          "refs"_a, "ref_to_genus"_a, "genusmat"_a)
      .def_prop_ro("ngenus", &dada2::TaxonomyModel::ngenus)
      .def_prop_ro("nlevel", &dada2::TaxonomyModel::nlevel);

  m.attr("NBOOT") = dada2::NBOOT;
  m.def(
      "assign_taxonomy",
      [](const dada2::TaxonomyModel &model, const std::vector<std::string> &seqs,
         const std::vector<std::string> &rcs, bool try_rc, uint64_t seed, int n_threads) {
        dada2::TaxonomyResult r;
        {
          PyHandlers handlers(std::nullopt);
          nb::gil_scoped_release release;
          r = dada2::assign_taxonomy(model, seqs, rcs, try_rc, seed, n_threads);
        }
        const size_t nseq = seqs.size();
        nb::dict d;
        d["tax"] = to_numpy(std::move(r.tax));
        d["boot"] = to_numpy(std::move(r.boot), {nseq, static_cast<size_t>(model.nlevel())});
        d["boot_tax"] = to_numpy(std::move(r.boot_tax), {nseq, static_cast<size_t>(dada2::NBOOT)});
        return d;
      },
      "model"_a, "seqs"_a, "rcs"_a, "try_rc"_a, "seed"_a, "n_threads"_a = 1,
      "Classify sequences; returns dict(tax, boot, boot_tax) with 0-based genus indices.");
  // ---------------------------------------------------------- ported kernels (spec §6) ---
  m.def("detect_quality_offset", &dada2::detect_quality_offset, "quals"_a,
        "ShortRead quality-type auto detection on a chunk: 64 (Solexa) or 33 (Sanger).");
  m.def(
      "filter_mate",
      [](const std::vector<std::string> &seqs, const std::vector<std::string> &quals, int trunc_q,
         int trunc_len, double max_len, int min_len, int trim_left, int trim_right, int max_n,
         int min_q, double max_ee, bool rm_phix, double rm_lowcomplex, int qual_offset,
         const std::string &phix, int n_threads) {
        dada2::MateFilterParams p;
        p.trunc_q = trunc_q; p.trunc_len = trunc_len; p.max_len = max_len; p.min_len = min_len;
        p.trim_left = trim_left; p.trim_right = trim_right; p.max_n = max_n; p.min_q = min_q;
        p.max_ee = max_ee; p.rm_phix = rm_phix; p.rm_lowcomplex = rm_lowcomplex;
        p.qual_offset = qual_offset;
        dada2::MateFilterResult r;
        {
          nb::gil_scoped_release release;
          r = dada2::filter_mate(seqs, quals, p, phix, n_threads);
        }
        std::vector<bool> pass(r.pass.begin(), r.pass.end());
        return nb::make_tuple(bools_to_numpy(pass), to_numpy(std::move(r.start)),
                              to_numpy(std::move(r.end)));
      },
      "seqs"_a, "quals"_a, "trunc_q"_a, "trunc_len"_a, "max_len"_a, "min_len"_a, "trim_left"_a,
      "trim_right"_a, "max_n"_a, "min_q"_a, "max_ee"_a, "rm_phix"_a, "rm_lowcomplex"_a,
      "qual_offset"_a, "phix"_a, "n_threads"_a = 1,
      "Per-mate filter/trim kernel (K1). Returns (pass, start, end) with 0-based [start, end).");
  m.def(
      "is_phix",
      [](const std::vector<std::string> &seqs, const std::string &phix, int word_size,
         int min_matches, bool non_overlapping) {
        std::vector<bool> r;
        {
          nb::gil_scoped_release release;
          r = dada2::is_phix(seqs, phix, word_size, min_matches, non_overlapping);
        }
        return bools_to_numpy(r);
      },
      "seqs"_a, "phix"_a, "word_size"_a = 16, "min_matches"_a = 2, "non_overlapping"_a = true,
      "Whether each sequence matches the phiX genome (or its reverse complement).");
  m.def(
      "seq_complexity",
      [](const std::vector<std::string> &seqs, int kmer_size, int window, int by, int n_threads) {
        std::vector<double> r;
        {
          nb::gil_scoped_release release;
          r = dada2::seq_complexity(seqs, kmer_size, window, by, n_threads);
        }
        return to_numpy(std::move(r));
      },
      "seqs"_a, "kmer_size"_a = 2, "window"_a = 0, "by"_a = 5, "n_threads"_a = 1,
      "Effective number of k-mers (exp Shannon entropy); window <= 0 uses the whole sequence (K6).");
  m.def("dna_less", &dada2::dna_less, "a"_a, "b"_a,
        "Biostrings lexical order of two DNA strings (IUPAC bit codes).");

  nb::class_<LockedDereplicator>(m, "Dereplicator", "Streaming dereplication (K2).")
      .def(nb::init<>())
      .def(
          "add_chunk",
          [](LockedDereplicator &d, const std::vector<std::string> &seqs,
             const std::vector<std::string> &quals, int qual_offset) {
            nb::gil_scoped_release release;
            std::lock_guard<std::mutex> lock(d.mutex);
            d.impl.add_chunk(seqs, quals, qual_offset);
          },
          "seqs"_a, "quals"_a, "qual_offset"_a)
      .def("finish", [](LockedDereplicator &d) {
        dada2::Dereplicator::Result r;
        {
          nb::gil_scoped_release release;
          std::lock_guard<std::mutex> lock(d.mutex);
          r = d.impl.finish();
        }
        const size_t nu = r.uniques.size();
        nb::dict out;
        out["uniques"] = nb::cast(std::move(r.uniques));
        out["abundance"] = to_numpy(std::move(r.abundance));
        out["quals"] = to_numpy(std::move(r.quals), {nu, static_cast<size_t>(r.maxlen)});
        out["map"] = to_numpy(std::move(r.map));
        out["n_zero_length"] = r.n_zero_length;
        return out;
      });

  m.def(
      "match_primer",
      [](const std::vector<std::string> &seqs, const std::string &primer, int max_mismatch,
         bool with_indels, bool fixed, int n_threads) {
        dada2::PrimerHits h;
        {
          nb::gil_scoped_release release;
          h = dada2::match_primer(seqs, primer, max_mismatch, with_indels, fixed, n_threads);
        }
        nb::dict d;
        d["nhits"] = to_numpy(std::move(h.nhits));
        d["first_start"] = to_numpy(std::move(h.first_start));
        d["first_end"] = to_numpy(std::move(h.first_end));
        d["last_start"] = to_numpy(std::move(h.last_start));
        d["last_end"] = to_numpy(std::move(h.last_end));
        return d;
      },
      "seqs"_a, "primer"_a, "max_mismatch"_a, "with_indels"_a, "fixed"_a, "n_threads"_a = 1,
      "Primer matches per sequence (K7): counts and first/last match (0-based start, end excl.).");
  m.def("match_pattern", &dada2::match_pattern, "subject"_a, "pattern"_a, "max_mismatch"_a,
        "with_indels"_a, "fixed"_a,
        "All matches (0-based start, exclusive end) of pattern in subject (Biostrings semantics).");
  m.def(
      "merge_pairs_core",
      [](const std::vector<std::string> &fwd, const std::vector<std::string> &rev_rc,
         const std::vector<int> &prefer, int match, int mismatch, int gap, bool trim_overhang,
         int n_threads) {
        dada2::MergeResult r;
        {
          nb::gil_scoped_release release;
          r = dada2::merge_pairs_core(fwd, rev_rc, prefer, match, mismatch, gap, trim_overhang,
                                      n_threads);
        }
        nb::dict d;
        d["nmatch"] = to_numpy(std::move(r.nmatch));
        d["nmismatch"] = to_numpy(std::move(r.nmismatch));
        d["nindel"] = to_numpy(std::move(r.nindel));
        d["sequence"] = nb::cast(std::move(r.sequence));
        return d;
      },
      "fwd"_a, "rev_rc"_a, "prefer"_a, "match"_a, "mismatch"_a, "gap"_a, "trim_overhang"_a,
      "n_threads"_a = 1, "Align, evaluate and merge read pairs (K3).");
  m.def(
      "collapse_no_mismatch",
      [](const std::vector<std::string> &seqs, int min_overlap, bool vec, int band, int n_threads) {
        std::vector<int> r;
        {
          nb::gil_scoped_release release;
          r = dada2::collapse_no_mismatch(seqs, min_overlap, vec, band, n_threads);
        }
        return to_numpy(std::move(r));
      },
      "seqs"_a, "min_overlap"_a, "vec"_a, "band"_a, "n_threads"_a = 1,
      "Greedy no-mismatch collapsing (K5): target index per sequence.");
  m.def(
      "is_shift_denovo",
      [](const std::vector<std::string> &seqs, const std::vector<int64_t> &abunds, int min_overlap,
         bool flag_subseqs, int n_threads) {
        std::vector<bool> r;
        {
          nb::gil_scoped_release release;
          r = dada2::is_shift_denovo(seqs, abunds, min_overlap, flag_subseqs, n_threads);
        }
        return bools_to_numpy(r);
      },
      "seqs"_a, "abunds"_a, "min_overlap"_a, "flag_subseqs"_a, "n_threads"_a = 1,
      "Shifted copies of more abundant sequences (K5).");
  m.def(
      "bimera_denovo",
      [](const std::vector<std::string> &seqs, const std::vector<double> &abunds, double min_fold,
         double min_parent_abund, bool allow_one_off, int min_one_off_par_dist, int match,
         int mismatch, int gap_p, int max_shift, int n_threads) {
        std::vector<bool> r;
        {
          nb::gil_scoped_release release;
          r = dada2::bimera_denovo(seqs, abunds, min_fold, min_parent_abund, allow_one_off,
                                   min_one_off_par_dist, match, mismatch, gap_p, max_shift,
                                   n_threads);
        }
        return bools_to_numpy(r);
      },
      "seqs"_a, "abunds"_a, "min_fold"_a, "min_parent_abund"_a, "allow_one_off"_a,
      "min_one_off_par_dist"_a, "match"_a, "mismatch"_a, "gap_p"_a, "max_shift"_a,
      "n_threads"_a = 1, "De novo bimera flags for a set of uniques (K4).");
  m.def(
      "species_hits",
      [](const std::vector<std::string> &queries, const std::vector<std::string> &refs, bool try_rc,
         int n_threads) {
        std::vector<std::vector<int>> r;
        {
          nb::gil_scoped_release release;
          r = dada2::species_hits(queries, refs, try_rc, n_threads);
        }
        return r;
      },
      "queries"_a, "refs"_a, "try_rc"_a, "n_threads"_a = 1,
      "Indices of the references containing each query exactly (K8).");
}
