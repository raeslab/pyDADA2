// K3: read-pair merging.
// Ported from DADA2 R/paired.R::mergePairs (the per-pair alignment, evaluation and consensus
// loop) (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// R aligns each unique (forward, rc(reverse)) pair with nwalign(band=-1, endsfree=TRUE) using
// MATCH=1, MISMATCH/GAP = -64 (maxMismatch == 0) or -8, then C_eval_pair and C_pair_consensus.
// pydada2: identical (forward, rc(reverse)) pairs, e.g. the same ASV pair in several samples of
// one call, are aligned and evaluated once; the consensus is formed per entry (prefer may differ).

#include <functional>
#include <string_view>
#include <unordered_map>

#include "../api.h"
#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

MergeResult merge_pairs_core(const std::vector<std::string> &fwd, const std::vector<std::string> &rev_rc,
                             const std::vector<int> &prefer, int match, int mismatch, int gap,
                             bool trim_overhang, int n_threads) {
  const size_t n = fwd.size();
  if (rev_rc.size() != n || prefer.size() != n)
    dada2_stop_invalid("Forward, reverse and prefer vectors must have the same length.");
  MergeResult r;
  r.nmatch.assign(n, 0);
  r.nmismatch.assign(n, 0);
  r.nindel.assign(n, 0);
  r.sequence.assign(n, std::string());

  // Distinct (forward, rc(reverse)) pairs, in order of first appearance.
  struct PairHash {
    size_t operator()(const std::pair<std::string_view, std::string_view> &p) const {
      const size_t h = std::hash<std::string_view>()(p.first);
      return h ^ (std::hash<std::string_view>()(p.second) + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2));
    }
  };
  std::unordered_map<std::pair<std::string_view, std::string_view>, size_t, PairHash> index;
  index.reserve(n);
  std::vector<size_t> pair_of(n), first;
  for (size_t i = 0; i < n; ++i) {
    auto ins = index.emplace(std::make_pair(std::string_view(fwd[i]), std::string_view(rev_rc[i])),
                             first.size());
    if (ins.second) first.push_back(i);
    pair_of[i] = ins.first->second;
  }
  std::vector<std::pair<std::string, std::string>> aligned(first.size());
  parallel_for(0, first.size(), [&](size_t b, size_t e) {
    for (size_t u = b; u < e; ++u) {
      const size_t i = first[u];
      aligned[u] = nwalign(fwd[i], rev_rc[i], match, mismatch, gap, gap, -1, true);
    }
  }, 8, n_threads);
  parallel_for(0, n, [&](size_t b, size_t e) {
    for (size_t i = b; i < e; ++i) {
      const auto &al = aligned[pair_of[i]];
      std::array<int, 3> ev = eval_pair(al.first, al.second);
      r.nmatch[i] = ev[0];
      r.nmismatch[i] = ev[1];
      r.nindel[i] = ev[2];
      r.sequence[i] = pair_consensus(al.first, al.second, prefer[i], trim_overhang);
    }
  }, 64, n_threads);
  return r;
}

}  // namespace dada2
