// K1: per-mate read filtering/trimming kernel.
// Ported from DADA2 R/filter.R::fastqFilter and ::fastqPairedFilter (and the ShortRead /
// Biostrings helpers they use: trimTails, narrow, .nFilter/alphabetFrequency) and
// R/filter.R::isPhiX (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// The R code applies each filter to the whole chunk in sequence; because every step only looks
// at the read itself, the same result is obtained by running the steps per read. For paired
// reads R keeps a pair only if both mates pass every step, i.e. pass_F && pass_R.
//
// The expected errors use a 256-entry table of the same pow(10, -q/10) values; the phiX check runs
// in parallel (match_ref).

#include <algorithm>
#include <cmath>
#include <unordered_set>

#include "../api.h"
#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

int detect_quality_offset(const std::vector<std::string> &quals) {
  // breaks: ':' (58) and 'J' (74) in the ASCII alphabet (.qualityTypeAuto)
  int mn = 256, mx = -1;
  const size_t n = std::min<size_t>(quals.size(), 10000);
  for (size_t i = 0; i < n; ++i) {
    for (unsigned char c : quals[i]) {
      if (c < mn) mn = c;
      if (c > mx) mx = c;
    }
  }
  if (mx >= 0 && mn > ':' && mx > 'J') return 64;
  return 33;
}

std::vector<bool> is_phix(const std::vector<std::string> &seqs, const std::string &phix,
                          int word_size, int min_matches, bool non_overlapping, int n_threads) {
  std::string rc(phix.rbegin(), phix.rend());
  for (char &c : rc) {
    switch (c) {
      case 'A': c = 'T'; break;
      case 'C': c = 'G'; break;
      case 'G': c = 'C'; break;
      case 'T': c = 'A'; break;
      default: break;
    }
  }
  std::vector<int> hits = match_ref(seqs, phix, word_size, non_overlapping, n_threads);
  std::vector<int> hits_rc = match_ref(seqs, rc, word_size, non_overlapping, n_threads);
  std::vector<bool> out(seqs.size());
  for (size_t i = 0; i < seqs.size(); ++i)
    out[i] = (hits[i] >= min_matches) || (hits_rc[i] >= min_matches);
  return out;
}

MateFilterResult filter_mate(const std::vector<std::string> &seqs,
                             const std::vector<std::string> &quals, const MateFilterParams &p,
                             const std::string &phix, int n_threads) {
  const size_t n = seqs.size();
  if (quals.size() != n) dada2_stop_invalid("Sequences and qualities differ in number.");
  MateFilterResult r;
  r.pass.assign(n, 0);
  r.start.assign(n, 0);
  r.end.assign(n, 0);

  // start/end as computed in R: start <- max(1, trimLeft+1); end <- truncLen;
  // if(end < start) end <- NA; end <- end - start + 1  (end relative to the trimmed read)
  const int start = std::max(1, p.trim_left + 1);
  int end_rel = -1;  // -1 == NA
  if (p.trunc_len >= start) end_rel = p.trunc_len - start + 1;
  const bool use_min_q = p.min_q > p.trunc_q;
  const bool use_ee = std::isfinite(p.max_ee);
  const unsigned char trunc_char = static_cast<unsigned char>(p.trunc_q + p.qual_offset);
  // Expected error of each quality character: the same pow() values, computed once.
  double ee_of_char[256];
  for (int c = 0; c < 256; ++c) {
    const int qq = c - p.qual_offset;
    ee_of_char[c] = pow(10.0, -qq / 10.0);
  }

  auto body = [&](size_t b, size_t e) {
    for (size_t i = b; i < e; ++i) {
      const std::string &s = seqs[i];
      const std::string &q = quals[i];
      if (q.size() != s.size()) dada2_stop_invalid("Sequence and quality lengths differ.");
      int len = static_cast<int>(s.size());
      // maxLen (raw reads)
      if (std::isfinite(p.max_len) && len > p.max_len) continue;
      // trimLeft: keep width >= start, then narrow(start)
      if (len < start) continue;
      int lo = start - 1;
      int hi = len;  // exclusive
      // trimRight: keep width > trimRight, then drop the last trimRight bases
      if (p.trim_right > 0) {
        if (hi - lo <= p.trim_right) continue;
        hi -= p.trim_right;
      }
      // truncQ: trimTails(k=1, a=truncQ): cut before the first base with qual char <= a;
      // zero-length results are removed
      for (int j = lo; j < hi; ++j) {
        if (static_cast<unsigned char>(q[j]) <= trunc_char) {
          hi = j;
          break;
        }
      }
      if (hi - lo <= 0) continue;
      // truncLen: keep width >= end, then narrow(1, end)
      if (end_rel >= 0) {
        if (hi - lo < end_rel) continue;
        hi = lo + end_rel;
      }
      // minLen
      if (hi - lo < p.min_len) continue;
      // maxN: number of non-ACGT letters
      int nother = 0;
      for (int j = lo; j < hi; ++j) {
        const char c = s[j];
        if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T')) ++nother;
      }
      if (nother > p.max_n) continue;
      // minQ: min(qual) > minQ (strict, as in R), only if minQ > truncQ
      if (use_min_q) {
        int mq = 1 << 30;
        for (int j = lo; j < hi; ++j) mq = std::min(mq, int(static_cast<unsigned char>(q[j])) - p.qual_offset);
        if (!(mq > p.min_q)) continue;
      }
      // maxEE (C_matrixEE)
      if (use_ee) {
        double ee = 0.0;
        for (int j = lo; j < hi; ++j) ee += ee_of_char[static_cast<unsigned char>(q[j])];
        if (!(ee <= p.max_ee)) continue;
      }
      r.start[i] = lo;
      r.end[i] = hi;
      r.pass[i] = 1;
    }
  };
  parallel_for(0, n, body, 1024, n_threads);

  // phiX and low complexity on the trimmed sequences of the reads that are still kept
  if (p.rm_phix || p.rm_lowcomplex != 0) {
    std::vector<size_t> idx;
    std::vector<std::string> trimmed;
    for (size_t i = 0; i < n; ++i) {
      if (r.pass[i]) {
        idx.push_back(i);
        trimmed.push_back(seqs[i].substr(r.start[i], r.end[i] - r.start[i]));
      }
    }
    if (p.rm_phix && !trimmed.empty()) {
      std::vector<bool> phi = is_phix(trimmed, phix, 16, 2, true, n_threads);
      for (size_t k = 0; k < idx.size(); ++k)
        if (phi[k]) r.pass[idx[k]] = 0;
    }
    if (p.rm_lowcomplex != 0 && !trimmed.empty()) {
      std::vector<double> cx = seq_complexity(trimmed, 2, 0, 5, n_threads);
      // R: is.lowc <- seqComplexity(...) < rm.lowcomplex (NA -> treated as low complexity)
      for (size_t k = 0; k < idx.size(); ++k)
        if (!(cx[k] >= p.rm_lowcomplex)) r.pass[idx[k]] = 0;
    }
  }
  return r;
}

}  // namespace dada2
