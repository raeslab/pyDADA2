// K8: exact-substring species matching.
// Ported from DADA2 R/taxonomy.R::assignSpecies (the Biostrings PDict/vcountPDict exact matching
// of each query against all reference sequences, and their reverse complements with tryRC)
// (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Instead of a PDict per batch of equal-length queries, each distinct query length L is handled
// with a hash table of the queries and one pass over all length-L windows of every reference.

#include <algorithm>
#include <unordered_map>

#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

namespace {

std::string rc_acgt(const std::string &s) {
  std::string r(s.rbegin(), s.rend());
  for (char &c : r) {
    switch (c) {
      case 'A': c = 'T'; break;
      case 'C': c = 'G'; break;
      case 'G': c = 'C'; break;
      case 'T': c = 'A'; break;
      case 'M': c = 'K'; break;
      case 'K': c = 'M'; break;
      case 'R': c = 'Y'; break;
      case 'Y': c = 'R'; break;
      case 'B': c = 'V'; break;
      case 'V': c = 'B'; break;
      case 'D': c = 'H'; break;
      case 'H': c = 'D'; break;
      default: break;  // S, W, N, '-' map to themselves
    }
  }
  return r;
}

}  // namespace

std::vector<std::vector<int>> species_hits(const std::vector<std::string> &queries,
                                           const std::vector<std::string> &refs, bool try_rc,
                                           int n_threads) {
  const size_t nq = queries.size();
  const size_t nr = refs.size();
  std::vector<std::string> rcs;
  if (try_rc) {
    rcs.reserve(nr);
    for (const auto &r : refs) rcs.push_back(rc_acgt(r));
  }
  // group queries by length
  std::unordered_map<size_t, std::unordered_map<std::string, std::vector<int>>> by_len;
  for (size_t i = 0; i < nq; ++i) by_len[queries[i].size()][queries[i]].push_back(static_cast<int>(i));
  std::vector<std::vector<int>> per_ref(nr);
  auto body = [&](size_t b, size_t e) {
    for (size_t r = b; r < e; ++r) {
      std::vector<int> &out = per_ref[r];
      for (const auto &kv : by_len) {
        const size_t L = kv.first;
        if (L == 0) continue;
        const auto &tab = kv.second;
        for (int pass = 0; pass < (try_rc ? 2 : 1); ++pass) {
          const std::string &s = pass == 0 ? refs[r] : rcs[r];
          if (s.size() < L) continue;
          for (size_t p = 0; p + L <= s.size(); ++p) {
            auto it = tab.find(s.substr(p, L));
            if (it != tab.end()) out.insert(out.end(), it->second.begin(), it->second.end());
          }
        }
      }
      std::sort(out.begin(), out.end());
      out.erase(std::unique(out.begin(), out.end()), out.end());
    }
  };
  parallel_for(0, nr, body, 16, n_threads);
  std::vector<std::vector<int>> hits(nq);
  for (size_t r = 0; r < nr; ++r)
    for (int q : per_ref[r]) hits[q].push_back(static_cast<int>(r));
  return hits;  // per query: sorted reference indices that contain it
}

}  // namespace dada2
