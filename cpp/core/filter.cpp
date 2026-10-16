// Adapted from DADA2 src/filter.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * C_matchRef -> dada2::match_ref, C_matrixEE -> dada2::matrix_ee (row-major input, NA_QUAL for NA).
//  * match_ref: 2-bit packed words with a bitmap prefilter for A/C/G/T references, in parallel.

#include <algorithm>
#include <unordered_set>

#include "api.h"
#include "dada.h"
#include "parallel.h"

namespace dada2 {

namespace {

int nt_code(char c) {
  switch(c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return -1;
  }
}

}  // namespace

// Counts word_size-mers of each sequence found in ref (circular), was C_matchRef.
// pydada2: when the reference is pure A/C/G/T and word_size <= 32 (always for isPhiX), words are
// compared as 2-bit codes (a 2^24-bit hashed bitmap rejects almost every window, a sorted array
// confirms hits) instead of an unordered_set<std::string> with a heap-allocated substr per
// window; a window containing another letter cannot equal a reference word. The scan, including
// the non_overlapping skip (j += word_size, then the loop's j++), is unchanged; sequences are
// processed in parallel.
std::vector<int> match_ref(const std::vector<std::string> &seqs, const std::string &ref_in,
                           int word_size_in, bool non_overlapping, int n_threads) {
  if(word_size_in < 1) { dada2_stop_invalid("word_size must be positive."); }
  const unsigned int word_size = word_size_in;
  std::vector<int> rval(seqs.size(), 0);

  bool packed = word_size <= 32;
  for(char c : ref_in) { if(nt_code(c) < 0) { packed = false; break; } }

  if(!packed) {  // upstream implementation
    std::string ref = ref_in;
    std::unordered_set<std::string> phash;
    unsigned int len = ref.size();
    ref.append(ref, 0, word_size);
    for(unsigned int i=0;i<len;i++) {
      phash.insert(ref.substr(i,word_size));
    }
    dada2::parallel_for(0, seqs.size(), [&](std::size_t begin, std::size_t end) {
      for(std::size_t i=begin;i<end;i++) {
        unsigned int slen = seqs[i].size();
        if(slen<word_size) { continue; }
        for(int j=0;j<=(int)(slen-word_size);j++) {
          if(phash.count(seqs[i].substr(j, word_size))) {
            rval[i]++;
            if(non_overlapping) { j+=word_size; }
          }
        }
      }
    }, 256, n_threads);
    return rval;
  }

  // Reference words (upstream: ref.substr(i, word_size) of the circularised reference, i < len;
  // a substring shorter than word_size, possible only when len < word_size, matches nothing).
  const uint64_t mask = word_size == 32 ? ~(uint64_t)0 : (((uint64_t)1 << (2*word_size)) - 1);
  std::string ref = ref_in;
  const size_t len = ref.size();
  ref.append(ref, 0, word_size);
  std::vector<uint64_t> words;
  words.reserve(len);
  for(size_t i=0;i<len && i+word_size<=ref.size();i++) {
    uint64_t v = 0;
    for(unsigned int q=0;q<word_size;q++) { v = (v << 2) | (uint64_t)nt_code(ref[i+q]); }
    words.push_back(v & mask);
  }
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  static constexpr unsigned int BITS = 24;  // static: MSVC rejects reading a local in a [] lambda
  std::vector<uint64_t> bitmap(((size_t)1 << BITS) / 64, 0);
  auto hash = [](uint64_t v) { return (uint32_t)((v * 0x9E3779B97F4A7C15ULL) >> (64 - BITS)); };
  for(uint64_t v : words) { const uint32_t h = hash(v); bitmap[h >> 6] |= (uint64_t)1 << (h & 63); }

  dada2::parallel_for(0, seqs.size(), [&](std::size_t begin, std::size_t end) {
    for(std::size_t i=begin;i<end;i++) {
      const std::string &sq = seqs[i];
      const size_t slen = sq.size();
      if(slen < word_size) { continue; }
      uint64_t v = 0;
      size_t run = 0;      // A/C/G/T letters ending at p
      size_t next_ok = 0;  // windows starting before this are skipped (non_overlapping)
      int count = 0;
      for(size_t p=0;p<slen;p++) {
        const int c = nt_code(sq[p]);
        if(c < 0) { run = 0; continue; }
        v = ((v << 2) | (uint64_t)c) & mask;
        if(++run < word_size) { continue; }
        const size_t j = p + 1 - word_size;  // window start
        if(j < next_ok) { continue; }
        const uint32_t h = hash(v);
        if(((bitmap[h >> 6] >> (h & 63)) & 1) && std::binary_search(words.begin(), words.end(), v)) {
          count++;
          if(non_overlapping) { next_ok = j + word_size + 1; }
        }
      }
      rval[i] = count;
    }
  }, 256, n_threads);
  return rval;
}

// Expected errors per row of a quality matrix; a NA_QUAL entry ends the row (was C_matrixEE).
std::vector<double> matrix_ee(const int32_t *inp, int nrow, int ncol) {
  int i,j;
  double ee;
  std::vector<double> rval(nrow);

  for(i=0;i<nrow;i++) {
    ee=0.0;
    for(j=0;j<ncol;j++) {
      int32_t q = inp[(size_t)i*ncol + j];
      if(q == NA_QUAL) { break; }
      ee += pow(10.0, -q/10.0);
    }
    rval[i] = ee;
  }
  return(rval);
}

}  // namespace dada2
