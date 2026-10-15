// K5: collapse sequences identical up to shifts/length, and shift detection.
// Ported from DADA2 R/multiSample.R::collapseNoMismatch and R/chimeras.R::isShiftDenovo,
// ::isShift, ::isShiftedPair, R/misc.R::nwhamming
// (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// collapse_no_mismatch precomputes the min_overlap prefixes once; its candidate blocks run on the
// shared thread pool (no thread creation per block).

#include <algorithm>

#include "../api.h"
#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

namespace {

constexpr int MATCH = 5, MISMATCH = -4, GAP = -8;  // getDadaOpt() defaults used by nwhamming

int nwhamming1(const std::string &a, const std::string &b, bool vec, int band) {
  std::pair<std::string, std::string> al;
  if (vec) {
    al = nwvec({a}, {b}, MATCH, MISMATCH, GAP, band, true)[0];
  } else {
    al = nwalign(a, b, MATCH, MISMATCH, GAP, GAP, band, true);
  }
  std::array<int, 3> ev = eval_pair(al.first, al.second);
  return ev[1] + ev[2];
}

}  // namespace

std::vector<int> collapse_no_mismatch(const std::vector<std::string> &seqs, int min_overlap,
                                      bool vec, int band, int n_threads) {
  // seqs are in R's processing order (decreasing total abundance). Returns, for each sequence,
  // the index (into seqs) of the output sequence it is merged into (itself if kept).
  const size_t n = seqs.size();
  std::vector<int> target(n, -1);
  std::vector<int> out;  // indices of kept sequences, in the order they were added
  // Prefixes of min_overlap characters (computed once, not once per pair).
  std::vector<std::string> prefixes(n);
  for (size_t q = 0; q < n; ++q) {
    prefixes[q] = seqs[q].substr(0, std::min<size_t>(seqs[q].size(), std::max(min_overlap, 0)));
  }
  const size_t block = static_cast<size_t>(std::max(1, n_threads)) * 4;
  std::vector<int> cand;
  for (size_t q = 0; q < n; ++q) {
    const std::string &query = seqs[q];
    const std::string &pq = prefixes[q];
    // Prescreen (cheap): refs sharing a prefix substring, in insertion order.
    cand.clear();
    for (int ref : out) {
      const std::string &rs = seqs[ref];
      if (rs.find(pq) != std::string::npos || query.find(prefixes[ref]) != std::string::npos) {
        cand.push_back(ref);
      }
    }
    int hit = -1;
    // Evaluate candidates in order, in parallel blocks; the first (lowest-position) hit wins.
    for (size_t b0 = 0; b0 < cand.size() && hit < 0; b0 += block) {
      const size_t b1 = std::min(cand.size(), b0 + block);
      std::vector<unsigned char> ok(b1 - b0, 0);
      auto body = [&](size_t b, size_t e) {
        for (size_t k = b; k < e; ++k) ok[k - b0] = nwhamming1(query, seqs[cand[k]], vec, band) == 0;
      };
      parallel_for(b0, b1, body, 1, n_threads);
      for (size_t k = b0; k < b1; ++k) {
        if (ok[k - b0]) {
          hit = cand[k];
          break;
        }
      }
    }
    if (hit >= 0) {
      target[q] = hit;
    } else {
      target[q] = static_cast<int>(q);
      out.push_back(static_cast<int>(q));
    }
  }
  return target;
}

std::vector<bool> is_shift_denovo(const std::vector<std::string> &seqs,
                                  const std::vector<int64_t> &abunds, int min_overlap,
                                  bool flag_subseqs, int n_threads) {
  const size_t n = seqs.size();
  if (abunds.size() != n) dada2_stop_invalid("Sequences and abundances differ in number.");
  std::vector<unsigned char> res(n, 0);
  auto body = [&](size_t b, size_t e) {
    for (size_t i = b; i < e; ++i) {
      const std::string &sq = seqs[i];
      for (size_t j = 0; j < n; ++j) {
        if (!(abunds[j] > abunds[i])) continue;  // parents: strictly more abundant
        const std::string &par = seqs[j];
        std::pair<std::string, std::string> al = nwalign(sq, par, MATCH, MISMATCH, GAP, GAP, -1, true);
        std::array<int, 3> ev = eval_pair(al.first, al.second);
        const int m = ev[0];
        if ((m < (int)sq.size() || flag_subseqs) && (m < (int)par.size() || flag_subseqs) &&
            m >= min_overlap && ev[1] == 0 && ev[2] == 0) {
          res[i] = 1;
          break;
        }
      }
    }
  };
  parallel_for(0, n, body, 1, n_threads);
  return std::vector<bool>(res.begin(), res.end());
}

}  // namespace dada2
