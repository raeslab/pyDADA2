// K3: read-pair merging.
// Ported from DADA2 R/paired.R::mergePairs (the per-pair alignment, evaluation and consensus
// loop) (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// R aligns each unique (forward, rc(reverse)) pair with nwalign(band=-1, endsfree=TRUE) using
// MATCH=1, MISMATCH/GAP = -64 (maxMismatch == 0) or -8, then C_eval_pair and C_pair_consensus.

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
  auto body = [&](size_t b, size_t e) {
    for (size_t i = b; i < e; ++i) {
      std::pair<std::string, std::string> al =
          nwalign(fwd[i], rev_rc[i], match, mismatch, gap, gap, -1, true);
      std::array<int, 3> ev = eval_pair(al.first, al.second);
      r.nmatch[i] = ev[0];
      r.nmismatch[i] = ev[1];
      r.nindel[i] = ev[2];
      r.sequence[i] = pair_consensus(al.first, al.second, prefer[i], trim_overhang);
    }
  };
  parallel_for(0, n, body, 8, n_threads);
  return r;
}

}  // namespace dada2
