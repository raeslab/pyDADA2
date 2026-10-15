// K8: exact-substring species matching.
// Ported from DADA2 R/taxonomy.R::assignSpecies (the Biostrings PDict/vcountPDict exact matching
// of each query against all reference sequences, and their reverse complements with tryRC)
// (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Instead of a PDict per batch of equal-length queries, each distinct query length L is handled
// with a hash table of the queries and one pass over all length-L windows of every reference.
// Windows are hashed with a rolling (Rabin-Karp) hash, so each pass is O(reference length);
// every hash hit is confirmed with memcmp, so the matches are exact.

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <unordered_map>

#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

namespace {

std::string rc_acgt(std::string_view s) {
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
                                           const std::vector<std::string_view> &refs,
                                           bool try_rc, int n_threads) {
  const size_t nq = queries.size();
  const size_t nr = refs.size();
  // Group the distinct queries by length: rolling hash of each -> (query text, query indices).
  constexpr uint64_t B = 1099511628211ULL;  // odd multiplier, arithmetic mod 2^64
  auto hash = [](const char *p, size_t n) {
    uint64_t h = 0;
    for (size_t i = 0; i < n; ++i) h = h * B + static_cast<unsigned char>(p[i]);
    return h;
  };
  struct Query {
    const std::string *seq;
    std::vector<int> idx;
  };
  struct LenGroup {
    size_t len;
    uint64_t pow;  // B^(len-1), to drop the leading character
    std::unordered_map<uint64_t, std::vector<Query>> tab;
    std::vector<uint64_t> filter;  // bitmap of hash bits, rejects most windows before tab.find
  };
  constexpr unsigned FILTER_BITS = 20;
  auto filter_slot = [](uint64_t h) { return h >> (64 - FILTER_BITS); };
  std::vector<LenGroup> groups;
  {
    std::unordered_map<size_t, size_t> gi;
    for (size_t i = 0; i < nq; ++i) {
      const std::string &q = queries[i];
      if (q.empty()) continue;
      auto ins = gi.emplace(q.size(), groups.size());
      if (ins.second) {
        uint64_t pw = 1;
        for (size_t k = 1; k < q.size(); ++k) pw *= B;
        groups.push_back({q.size(), pw, {}, std::vector<uint64_t>(1u << (FILTER_BITS - 6))});
      }
      LenGroup &g = groups[ins.first->second];
      const uint64_t h = hash(q.data(), q.size());
      g.filter[filter_slot(h) >> 6] |= uint64_t{1} << (filter_slot(h) & 63);
      auto &bucket = g.tab[h];
      auto it = std::find_if(bucket.begin(), bucket.end(),
                             [&](const Query &x) { return *x.seq == q; });
      if (it == bucket.end()) bucket.push_back({&q, {static_cast<int>(i)}});
      else it->idx.push_back(static_cast<int>(i));
    }
  }
  std::vector<std::vector<int>> per_ref(nr);
  auto body = [&](size_t b, size_t e) {
    std::string rc;  // reverse complement of the current reference (try_rc), built per thread
    for (size_t r = b; r < e; ++r) {
      std::vector<int> &out = per_ref[r];
      if (try_rc) rc = rc_acgt(refs[r]);
      for (const auto &g : groups) {
        const size_t L = g.len;
        for (int pass = 0; pass < (try_rc ? 2 : 1); ++pass) {
          const std::string_view s = pass == 0 ? refs[r] : std::string_view(rc);
          if (s.size() < L) continue;
          const char *d = s.data();
          uint64_t h = hash(d, L);
          for (size_t p = 0;; ++p) {
            const uint64_t slot = filter_slot(h);
            if (g.filter[slot >> 6] >> (slot & 63) & 1) {
              auto it = g.tab.find(h);
              if (it != g.tab.end())
                for (const auto &q : it->second)
                  if (std::memcmp(q.seq->data(), d + p, L) == 0)
                    out.insert(out.end(), q.idx.begin(), q.idx.end());
            }
            if (p + L >= s.size()) break;
            h = (h - static_cast<unsigned char>(d[p]) * g.pow) * B +
                static_cast<unsigned char>(d[p + L]);
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
