// K6: sequence complexity.
// Ported from DADA2 R/filter.R::seqComplexity and ::sindex (with Biostrings
// oligonucleotideFrequency semantics) (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// The entropy terms are summed like R's sum(), which accumulates doubles in a long double
// (LDOUBLE) and rounds the total to double once: with plain double accumulation a uniform
// spectrum such as NGAGTACCNN (6 distinct 2-mers) gives 5.999999999999998 instead of R's 6, which
// flips rm_lowcomplex thresholds. Where long double has no extended precision (MSVC, Apple
// arm64) a double-double (TwoSum) accumulator stands in for R's 80-bit sum (the fixtures come
// from R on x86_64 Linux).

#include <cfloat>
#include <cmath>

#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

namespace {

inline int nt_index(char c) {
  switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return -1;
  }
}

// Accumulator reproducing R's sum() of doubles (long double accumulation, rounded once).
#if LDBL_MANT_DIG < 64
struct ExtendedSum {  // long double == double (MSVC, Apple arm64): double-double (TwoSum) instead
  double hi = 0.0, lo = 0.0;
  void add(double x) {
    const double s = hi + x;
    const double bp = s - hi;
    lo += (hi - (s - bp)) + (x - bp);
    hi = s;
  }
  double value() const { return hi + lo; }
};
#else
struct ExtendedSum {
  long double acc = 0.0L;
  void add(double x) { acc += x; }
  double value() const { return static_cast<double>(acc); }
};
#endif

// exp(Shannon entropy) of the k-mer counts of s[lo, hi); k-mers with non-ACGT letters are not
// counted (oligonucleotideFrequency). Summation follows R: k-mers in lexicographic order.
double sindex(const std::string &s, size_t lo, size_t hi, int k, std::vector<double> &counts) {
  std::fill(counts.begin(), counts.end(), 0.0);
  double total = 0.0;
  if (hi >= lo + static_cast<size_t>(k)) {
    for (size_t i = lo; i + k <= hi; ++i) {
      int code = 0;
      bool ok = true;
      for (int j = 0; j < k; ++j) {
        const int v = nt_index(s[i + j]);
        if (v < 0) {
          ok = false;
          break;
        }
        code = code * 4 + v;
      }
      if (ok) {
        counts[code] += 1.0;
        total += 1.0;
      }
    }
  }
  if (total == 0.0) return std::nan("");  // R: 0/0 -> NaN -> NA
  ExtendedSum acc;
  for (double c : counts) {
    if (c > 0) {
      const double y = c / total;
      acc.add(-y * std::log(y));
    }
  }
  return std::exp(acc.value());
}

}  // namespace

std::vector<double> seq_complexity(const std::vector<std::string> &seqs, int kmer_size, int window,
                                   int by, int n_threads) {
  if (kmer_size < 1 || kmer_size > 12) dada2_stop_invalid("kmer_size must be between 1 and 12.");
  if (window > 0 && kmer_size >= window)
    dada2_stop_invalid("The window over which kmer frequency is calculated must be larger than the kmerSize.");
  if (window > 0 && by < 1) dada2_stop_invalid("by must be positive.");
  const size_t nk = static_cast<size_t>(1) << (2 * kmer_size);
  std::vector<double> out(seqs.size());
  size_t maxw = 0;
  for (const auto &s : seqs) maxw = std::max(maxw, s.size());

  auto body = [&](size_t b, size_t e) {
    std::vector<double> counts(nk);
    for (size_t i = b; i < e; ++i) {
      const std::string &s = seqs[i];
      if (window <= 0) {
        out[i] = sindex(s, 0, s.size(), kmer_size, counts);
        continue;
      }
      // R: si <- rep(4^k); for(i in seq(1, max(width)-window, by=by)) {
      //      keep <- width >= i+window-1; si[keep] <- pmin(si[keep], sindex(window)) }
      double si = static_cast<double>(nk);
      for (long long st = 1; st <= static_cast<long long>(maxw) - window; st += by) {
        if (static_cast<long long>(s.size()) >= st + window - 1) {
          const double v = sindex(s, static_cast<size_t>(st - 1), static_cast<size_t>(st - 1 + window),
                                  kmer_size, counts);
          if (std::isnan(v) || std::isnan(si)) si = std::nan("");  // pmin propagates NA
          else if (v < si) si = v;
        }
      }
      out[i] = si;
    }
  };
  parallel_for(0, seqs.size(), body, 256, n_threads);
  return out;
}

}  // namespace dada2
