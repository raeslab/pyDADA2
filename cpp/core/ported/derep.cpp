// K2: streaming dereplication.
// Ported from DADA2 R/sequenceIO.R::derepFastq and ::qtables2 (and ShortRead srsort/srrank)
// (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// R semantics reproduced:
//  * per chunk: zero-length reads are dropped (map = NA); the chunk's distinct sequences are taken
//    in ShortRead lexical order (srsort), new ones are appended to the running uniques in that order;
//  * qualities are summed per unique (rowsum) and averaged at the end;
//  * final order: stable sort by decreasing abundance (R's order()).

#include <algorithm>
#include <cmath>
#include <numeric>

#include "../dada.h"
#include "ported.h"

namespace dada2 {

namespace {

// Biostrings DNA_CODES: one bit per base; ambiguity codes are unions. Unknown letters sort last.
inline int dna_code(unsigned char c) {
  switch (c) {
    case 'A': return 1;
    case 'C': return 2;
    case 'M': return 3;
    case 'G': return 4;
    case 'R': return 5;
    case 'S': return 6;
    case 'V': return 7;
    case 'T': return 8;
    case 'W': return 9;
    case 'Y': return 10;
    case 'H': return 11;
    case 'K': return 12;
    case 'D': return 13;
    case 'B': return 14;
    case 'N': return 15;
    case '-': return 16;
    case '+': return 32;
    case '.': return 64;
    default: return 128 + c;
  }
}

}  // namespace

bool dna_less(const std::string &a, const std::string &b) {
  const size_t n = std::min(a.size(), b.size());
  for (size_t i = 0; i < n; ++i) {
    const int ca = dna_code(static_cast<unsigned char>(a[i]));
    const int cb = dna_code(static_cast<unsigned char>(b[i]));
    if (ca != cb) return ca < cb;
  }
  return a.size() < b.size();
}

void Dereplicator::add_chunk(const std::vector<std::string> &seqs,
                             const std::vector<std::string> &quals, int qual_offset) {
  const size_t n = seqs.size();
  if (quals.size() != n) dada2_stop_invalid("Sequences and qualities differ in number.");
  // Distinct non-empty sequences of this chunk, in lexical order.
  std::vector<const std::string *> distinct;
  {
    std::unordered_map<std::string, int> seen;
    seen.reserve(n * 2);
    for (size_t i = 0; i < n; ++i) {
      if (seqs[i].empty()) continue;
      if (seen.emplace(seqs[i], 0).second) distinct.push_back(&seqs[i]);
    }
    if (distinct.empty() && n > 0) {
      dada2_stop_invalid("Only zero-length sequences detected during dereplication.");
    }
  }
  std::sort(distinct.begin(), distinct.end(),
            [](const std::string *a, const std::string *b) { return dna_less(*a, *b); });
  for (const std::string *s : distinct) {
    if (index_.find(*s) == index_.end()) {
      index_.emplace(*s, static_cast<int64_t>(uniques_.size()));
      uniques_.push_back(*s);
      counts_.push_back(0);
      qsums_.emplace_back(s->size(), 0.0);
      maxlen_ = std::max(maxlen_, static_cast<int>(s->size()));
    }
  }
  for (size_t i = 0; i < n; ++i) {
    const std::string &s = seqs[i];
    if (s.empty()) {
      map_.push_back(-1);
      ++n_zero_;
      continue;
    }
    if (quals[i].size() != s.size()) dada2_stop_invalid("Sequence and quality lengths differ.");
    const int64_t u = index_.find(s)->second;
    map_.push_back(u);
    counts_[u] += 1;
    std::vector<double> &qs = qsums_[u];
    const std::string &q = quals[i];
    for (size_t j = 0; j < q.size(); ++j) {
      qs[j] += static_cast<double>(static_cast<int>(static_cast<unsigned char>(q[j])) - qual_offset);
    }
  }
}

Dereplicator::Result Dereplicator::finish() const {
  Result r;
  const size_t nu = uniques_.size();
  std::vector<size_t> ord(nu);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(),
                   [&](size_t a, size_t b) { return counts_[a] > counts_[b]; });
  std::vector<int64_t> rank(nu);
  for (size_t k = 0; k < nu; ++k) rank[ord[k]] = static_cast<int64_t>(k);

  r.maxlen = maxlen_;
  r.uniques.reserve(nu);
  r.abundance.reserve(nu);
  r.quals.assign(nu * static_cast<size_t>(maxlen_), std::nan(""));
  for (size_t k = 0; k < nu; ++k) {
    const size_t u = ord[k];
    r.uniques.push_back(uniques_[u]);
    r.abundance.push_back(counts_[u]);
    const double c = static_cast<double>(counts_[u]);
    const std::vector<double> &qs = qsums_[u];
    for (size_t j = 0; j < qs.size(); ++j) r.quals[k * maxlen_ + j] = qs[j] / c;
  }
  r.map.reserve(map_.size());
  for (int64_t m : map_) r.map.push_back(m < 0 ? -1 : rank[m]);
  r.n_zero_length = n_zero_;
  return r;
}

}  // namespace dada2
