// K4: de novo bimera detection over a set of unique sequences.
// Ported from DADA2 R/chimeras.R::isBimeraDenovo (the per-sequence parent selection and
// isBimera loop) (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// The parents are passed to is_bimera by pointer (they were copied for every query sequence).

#include "../api.h"
#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

std::vector<bool> bimera_denovo(const std::vector<std::string> &seqs,
                                const std::vector<double> &abunds, double min_fold,
                                double min_parent_abund, bool allow_one_off,
                                int min_one_off_par_dist, int match, int mismatch, int gap_p,
                                int max_shift, int n_threads) {
  const size_t n = seqs.size();
  if (abunds.size() != n) dada2_stop_invalid("Sequences and abundances differ in number.");
  std::vector<unsigned char> res(n, 0);
  auto body = [&](size_t b, size_t e) {
    std::vector<const std::string *> pars;
    for (size_t i = b; i < e; ++i) {
      // R: pars <- names(unqs)[(unqs > minFold*abund & unqs > minParentAbundance)]
      pars.clear();
      const double a = abunds[i];
      for (size_t j = 0; j < n; ++j) {
        if (abunds[j] > min_fold * a && abunds[j] > min_parent_abund) pars.push_back(&seqs[j]);
      }
      if (pars.size() < 2) continue;
      res[i] = is_bimera_ptr(seqs[i], pars, allow_one_off, min_one_off_par_dist, match, mismatch,
                             gap_p, max_shift);
    }
  };
  parallel_for(0, n, body, 1, n_threads);
  return std::vector<bool>(res.begin(), res.end());
}

}  // namespace dada2
