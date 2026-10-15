// K7: primer matching.
// Ported from DADA2 R/filter.R::removePrimers, reproducing the semantics of the Biostrings
// matchPattern/vmatchPattern calls it makes (https://github.com/benjjneb/dada2, commit b61fe2e,
// v1.36.0). Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Mismatch mode (with.indels=FALSE): every alignment start s in [-k, n-m+k] with at most k
// mismatches is a match; pattern positions outside the subject count as mismatches
// ("out-of-limits" matches, as reported by Biostrings).
//
// Indel mode (with.indels=TRUE): Biostrings documents "best local matches", but its implementation
// (src/match_pattern_indels.c, lowlevel_matching.c; Artistic-2.0) differs from that definition in
// edge cases (large max.mismatch, subject ends). The algorithm used by Biostrings is
// re-implemented here (not copied) so that results are identical to R's removePrimers.

#include <algorithm>

#include "../dada.h"
#include "../parallel.h"
#include "ported.h"

namespace dada2 {

namespace {

inline int iupac_bits(unsigned char c) {
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
    default: return 0;
  }
}

struct Matcher {
  const std::string &pat;
  bool fixed;
  std::vector<int> pbits;
  Matcher(const std::string &p, bool f) : pat(p), fixed(f), pbits(p.size()) {
    for (size_t i = 0; i < p.size(); ++i) pbits[i] = iupac_bits(static_cast<unsigned char>(p[i]));
  }
  bool eq(size_t j, char c) const {
    if (fixed) return pat[j] == c;
    return (pbits[j] & iupac_bits(static_cast<unsigned char>(c))) != 0;
  }
};

void mismatch_matches(const std::string &s, const Matcher &mt, int k,
                      std::vector<std::pair<int, int>> &out) {
  const int n = static_cast<int>(s.size());
  const int m = static_cast<int>(mt.pat.size());
  for (int st = -k; st <= n - m + k; ++st) {
    int mm = 0;
    for (int j = 0; j < m && mm <= k; ++j) {
      const int pos = st + j;
      if (pos < 0 || pos >= n || !mt.eq(j, s[pos])) ++mm;
    }
    if (mm <= k) out.emplace_back(st, st + m);
  }
}

// Banded edit distance of P1 against S starting at Ploffset, with a free end in S (Biostrings
// lowlevel _nedit_for_Ploffset). Returns the minimal edit distance and sets min_width to the
// shortest width of S achieving it. Positions beyond S count as mismatches. P1 is the pattern
// suffix pat[p1_off, m). min_width is left unchanged when P1 is empty (as upstream).
int nedit_for_ploffset(const Matcher &mt, int p1_off, const std::string &s, int ploffset,
                       int max_nedit, int &min_width) {
  const int m = static_cast<int>(mt.pat.size());
  const int pl = m - p1_off;
  const int n = static_cast<int>(s.size());
  if (pl == 0) return 0;
  const int max_nedit_plus1 = max_nedit + 1;
  if (max_nedit > pl) max_nedit = pl;
  const int row_length = 2 * max_nedit + 1;
  std::vector<int> buf1(row_length, 0), buf2(row_length, 0);
  int *prev_row = buf1.data(), *curr_row = buf2.data();
  int min_si = ploffset;
  int a, B, b, pi, si, min_nedit;
  auto propagate = [&](int B, int si, int pi) {
    const bool ok = si >= 0 && si < n && mt.eq(p1_off + pi, s[si]);
    int nedit = prev_row[B] + (ok ? 0 : 1);
    int B2;
    if ((B2 = B - 1) >= 0 && curr_row[B2] + 1 < nedit) nedit = curr_row[B2] + 1;
    if ((B2 = B + 1) < row_length && prev_row[B2] + 1 < nedit) nedit = prev_row[B2] + 1;
    curr_row[B] = nedit;
  };
  // stage 0
  for (B = max_nedit, b = 0; B < row_length; B++, b++) curr_row[B] = b;
  // stage 1
  for (a = 1, pi = 0; a < max_nedit; a++, pi++) {
    std::swap(prev_row, curr_row);
    B = max_nedit - a;
    curr_row[B++] = a;
    for (si = min_si; B < row_length; B++, si++) propagate(B, si, pi);
  }
  // stage 2
  std::swap(prev_row, curr_row);
  B = 0;
  curr_row[B++] = min_nedit = a;
  min_width = 0;
  for (si = min_si; B < row_length; B++, si++) {
    propagate(B, si, pi);
    if (curr_row[B] < min_nedit) {
      min_nedit = curr_row[B];
      min_width = si - ploffset + 1;
    }
  }
  a++;
  pi++;
  // stage 3 (with bailout)
  for (; pi < pl; pi++, a++, min_si++) {
    std::swap(prev_row, curr_row);
    min_nedit = a;
    min_width = 0;
    for (B = 0, si = min_si; B < row_length; B++, si++) {
      propagate(B, si, pi);
      if (curr_row[B] < min_nedit) {
        min_nedit = curr_row[B];
        min_width = si - ploffset + 1;
      }
    }
    if (min_nedit >= max_nedit_plus1) break;
  }
  return min_nedit;
}

// Mismatches of pat[p1_off, m) aligned at S[pshift...] with early stop after max_nmis
// (Biostrings _nmismatch_at_Pshift).
int nmismatch_at_pshift(const Matcher &mt, int p1_off, const std::string &s, int pshift,
                        int max_nmis) {
  const int m = static_cast<int>(mt.pat.size());
  const int n = static_cast<int>(s.size());
  int nmis = 0;
  for (int i = p1_off, j = pshift; i < m; ++i, ++j) {
    if (j >= 0 && j < n && mt.eq(i, s[j])) continue;
    if (nmis++ >= max_nmis) break;
  }
  return nmis;
}

// Biostrings matchPattern(with.indels=TRUE): scans S; each position whose letter matches some
// pattern letter anchors the first such pattern offset i0 (the pattern prefix before i0 counts
// as i0 edits); the rest of the pattern is aligned with a free end. Overlapping candidates are
// resolved with Biostrings' "provisory match" rule.
void indel_matches(const std::string &s, const Matcher &mt, int k,
                   std::vector<std::pair<int, int>> &out) {
  const int n = static_cast<int>(s.size());
  const int m = static_cast<int>(mt.pat.size());
  if (m == 0) return;
  // byte -> first pattern offset matching it
  int byte2offset[256];
  for (int y = 0; y < 256; ++y) {
    byte2offset[y] = -1;
    for (int i = 0; i < m; ++i) {
      if (mt.eq(i, static_cast<char>(y))) {
        byte2offset[y] = i;
        break;
      }
    }
  }
  bool have_prov = false;
  int prov_start = 0, prov_end = 0, prov_width = 0, prov_nedit = 0;  // 0-based start
  auto report = [&](int start, int width) { out.emplace_back(start, start + width); };
  auto provisory = [&](int start, int width, int nedit) {
    const int end = start + width - 1;
    if (have_prov) {
      if (end > prov_end) report(prov_start, prov_width);
      else if (nedit > prov_nedit) return;
    }
    have_prov = true;
    prov_start = start;
    prov_end = end;
    prov_width = width;
    prov_nedit = nedit;
  };
  int width1 = 0;  // as upstream, keeps its previous value when the pattern suffix is empty
  for (int j0 = 0; j0 < n; ++j0) {
    const int i0 = byte2offset[static_cast<unsigned char>(s[j0])];
    if (i0 < 0) continue;
    const int max_nmis1 = k - i0;
    if (max_nmis1 < 0) continue;
    int nedit1;
    if (max_nmis1 == 0) {
      nedit1 = nmismatch_at_pshift(mt, i0 + 1, s, j0 + 1, max_nmis1);
      width1 = m - i0 - 1;
    } else {
      nedit1 = nedit_for_ploffset(mt, i0 + 1, s, j0 + 1, max_nmis1, width1);
    }
    if (nedit1 <= max_nmis1) provisory(j0, width1 + 1, nedit1 + i0);
  }
  if (have_prov) report(prov_start, prov_width);
}

// Biostrings with max.mismatch >= pattern length (both modes): every window of the pattern's
// width that overlaps the subject, i.e. starts -(m-1) .. n-1 (0-based), is reported. Verified
// against Biostrings on 662 random cases (tests/reference/data/match_pattern_degenerate.json).
bool all_windows(int n, int m, int k, std::vector<std::pair<int, int>> &out) {
  if (m < 1 || k < m) return false;
  for (int st = 1 - m; st <= n - 1; ++st) out.emplace_back(st, st + m);
  return true;
}

}  // namespace

std::vector<std::pair<int, int>> match_pattern(const std::string &subject,
                                               const std::string &pattern, int max_mismatch,
                                               bool with_indels, bool fixed) {
  if (max_mismatch < 0) dada2_stop_invalid("max_mismatch must be >= 0.");
  Matcher mt(pattern, fixed);
  std::vector<std::pair<int, int>> out;
  if (all_windows(static_cast<int>(subject.size()), static_cast<int>(pattern.size()), max_mismatch, out))
    return out;
  if (with_indels) indel_matches(subject, mt, max_mismatch, out);
  else mismatch_matches(subject, mt, max_mismatch, out);
  return out;
}

PrimerHits match_primer(const std::vector<std::string> &seqs, const std::string &primer,
                        int max_mismatch, bool with_indels, bool fixed, int n_threads) {
  if (primer.empty()) dada2_stop_invalid("Primer must not be empty.");
  if (max_mismatch < 0) dada2_stop_invalid("max_mismatch must be >= 0.");
  const size_t n = seqs.size();
  PrimerHits h;
  h.nhits.assign(n, 0);
  h.first_start.assign(n, 0);
  h.first_end.assign(n, 0);
  h.last_start.assign(n, 0);
  h.last_end.assign(n, 0);
  Matcher mt(primer, fixed);
  auto body = [&](size_t b, size_t e) {
    std::vector<std::pair<int, int>> hits;
    for (size_t i = b; i < e; ++i) {
      hits.clear();
      if (all_windows(static_cast<int>(seqs[i].size()), static_cast<int>(primer.size()), max_mismatch, hits)) {
      } else if (with_indels) {
        indel_matches(seqs[i], mt, max_mismatch, hits);
      } else {
        mismatch_matches(seqs[i], mt, max_mismatch, hits);
      }
      h.nhits[i] = static_cast<int>(hits.size());
      if (!hits.empty()) {
        h.first_start[i] = hits.front().first;
        h.first_end[i] = hits.front().second;
        h.last_start[i] = hits.back().first;
        h.last_end[i] = hits.back().second;
      }
    }
  };
  parallel_for(0, n, body, 64, n_threads);
  return h;
}

}  // namespace dada2
