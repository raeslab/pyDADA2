// Adapted from DADA2 src/evaluate.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * Rcpp exports -> dada2:: API functions with plain C++ types (api.h).
//  * Unequal-length alignment strings: upstream warned and returned NULL; now
//    std::invalid_argument. Empty strings in eval_pair return zeros (upstream read out of bounds).
//  * Temporary C buffers are managed by std::vector/std::unique_ptr (exception safety).
//  * nwalign rejects non-ACGT input (R checks this in its nwalign() wrapper; the C code indexes
//    the 4x4 score matrix with the nucleotide code and read out of bounds otherwise).
//  * k-mer helpers validate kmer_size before sizing their 4^k buffers (upstream computed
//    1 << 2k first: undefined for k >= 16, a 2 x 2 GiB allocation for k = 15).
//  * kord_dist(SSE=1) uses kord_dist_SSEi on every platform (exact scalar emulation off x86).

#include <memory>

#include "api.h"
#include "dada.h"

namespace {

struct CFree {
  void operator()(void *p) const { free(p); }
};
using CStr = std::unique_ptr<char, CFree>;

// The k-mer routines (kmers.cpp) accept 3..8 (count vectors) or 1..8 (ordered k-mers).
void check_kmer_size(int kmer_size, int min_k) {
  if(kmer_size < min_k || kmer_size > 8) { dada2_stop_invalid("Invalid kmer-size (must be %i..8).", min_k); }
}

struct AlFree {
  void operator()(char **al) const { al_free(al); }
};
using Al = std::unique_ptr<char *, AlFree>;

}  // namespace

namespace dada2 {

//------------------------------------------------------------------
// Exposes Needleman-Wunsch alignment (was C_nwalign).
std::pair<std::string, std::string> nwalign(const std::string &s1, const std::string &s2, int match,
                                            int mismatch, int gap_p, int homo_gap_p, int band,
                                            bool endsfree) {
  int i, j;
  for(bool ok : is_acgt({s1, s2})) {
    if(!ok) { dada2_stop_invalid("Sequences must contain only A/C/G/T characters."); }
  }
  // Make integer-ized c-style sequence strings
  CStr seq1((char *) malloc(s1.length()+1)); //E
  CStr seq2((char *) malloc(s2.length()+1)); //E
  if (!seq1 || !seq2)  dada2_stop("Memory allocation failed.");
  nt2int(seq1.get(), s1.c_str());
  nt2int(seq2.get(), s2.c_str());
  // Make  c-style 2d score array
  int c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      if(i==j) { c_score[i][j] = match; }
      else { c_score[i][j] = mismatch; }
    }
  }
  // Perform alignment and convert back to ACGT
  Al al;
  if(endsfree) {
    if(gap_p == homo_gap_p) {
      al.reset(nwalign_endsfree(seq1.get(), s1.length(), seq2.get(), s2.length(), c_score, gap_p, band));
    } else {
      al.reset(nwalign_endsfree_homo(seq1.get(), s1.length(), seq2.get(), s2.length(), c_score, gap_p, homo_gap_p, band));
    }
  } else {
    if(gap_p != homo_gap_p) {
      dada2_log("Warning: A separate homopolymer gap penalty isn't implemented when endsfree=FALSE.\n\tAll gaps will be penalized by the regular gap penalty.\n");
    }
    al.reset(::nwalign(seq1.get(), s1.length(), seq2.get(), s2.length(), c_score, gap_p, band));
  }
  int2nt(al.get()[0], al.get()[0]);
  int2nt(al.get()[1], al.get()[1]);
  return {std::string(al.get()[0]), std::string(al.get()[1])};
}

//------------------------------------------------------------------
// Calculates the number of matches/mismatches/internal_indels in an alignment (was C_eval_pair).
// Ignores indels at the ends of the alignment.
std::array<int, 3> eval_pair(const std::string &s1, const std::string &s2) {
  int match, mismatch, indel, start, end;
  bool s1gap, s2gap;
  if(s1.size() != s2.size()) {
    dada2_stop_invalid("Aligned strings are not the same length.");
  }
  if(s1.empty()) { return {0, 0, 0}; }

  // Find start of align (end of initial gapping)
  s1gap = s2gap = true;
  start = -1;
  do {
    start++;
    s1gap = s1gap && (s1[start] == '-');
    s2gap = s2gap && (s2[start] == '-');
  } while((s1gap || s2gap) && start<(int)s1.size());

  // Find end of align (start of terminal gapping)
  s1gap = s2gap = true;
  end = s1.size();
  do {
    end--;
    s1gap = s1gap && (s1[end] == '-');
    s2gap = s2gap && (s2[end] == '-');
  } while((s1gap || s2gap) && end>=start);

  // Count the matches, mismatches, indels within the internal part of alignment.
  match = mismatch = indel = 0;
  for(int i=start;i<=end;i++) {
    if(s1[i]=='-' || s2[i]=='-') {
      indel++;
    } else if(s1[i] == s2[i]) {
      match++;
    } else {
      mismatch++;
    }
  }
  return {match, mismatch, indel};
}

//------------------------------------------------------------------
// Calculates the consensus of two sequences (prefer sequence wins mismatches) (was C_pair_consensus).
std::string pair_consensus(const std::string &s1, const std::string &s2, int prefer,
                           bool trim_overhang) {
  int i;
  if(s1.size() != s2.size()) {
    dada2_stop_invalid("Aligned strings are not the same length.");
  }

  std::string oseq(s1.size(), '\0');
  for(i=0;i<(int)s1.size();i++) {
    if(s1[i] == s2[i]) {
      oseq[i] = s1[i];
    } else if(s2[i] == '-') {
      oseq[i] = s1[i];
    } else if(s1[i] == '-') {
      oseq[i] = s2[i];
    } else {
      if(prefer==1) {
        oseq[i] = s1[i]; // s1 wins mismatches
      } else if(prefer==2) {
        oseq[i] = s2[i];
      } else {
        oseq[i] = 'N'; // Should never happen
      }
    }
  }
  // Trim overhangs
  if(trim_overhang) {
    for(i=0;i<(int)s1.size();i++) {
      if(s1[i] != '-') { break; }
      oseq[i] = '-';
    }
    for(i=s1.size()-1;i>=0;i--) {
      if(s2[i] != '-') { break; }
      oseq[i] = '-';
    }
  }
  // Remove remaining gaps
  std::string out;
  out.reserve(oseq.size());
  for(i=0;i<(int)s1.size();i++) {
    if(oseq[i] != '-') {
      out.push_back(oseq[i]);
    }
  }
  return out;
}

//------------------------------------------------------------------
// Checks a vector of character sequences for whether they are entirely ACGT (was C_isACGT).
std::vector<bool> is_acgt(const std::vector<std::string> &seqs) {
  std::vector<bool> isACGT(seqs.size());
  for(size_t i=0; i<seqs.size(); i++) {
    bool justACGT = true;
    for(char c : seqs[i]) {
      if(!(c == 'A' || c == 'C' || c == 'G' || c == 'T')) {
        justACGT = false;
        break;
      }
    }
    isACGT[i] = justACGT;
  }
  return isACGT;
}

std::vector<double> kmer_dist(const std::vector<std::string> &s1, const std::vector<std::string> &s2,
                              int kmer_size) {
  size_t len1 = 0, len2 = 0;
  check_kmer_size(kmer_size, 3);
  size_t n_kmers = (1 << (2*kmer_size));  // 4^k kmers

  size_t nseqs = s1.size();
  if(nseqs != s2.size()) { dada2_stop_invalid("Mismatched numbers of sequences."); }

  std::vector<double> kdist(nseqs);
  std::vector<uint16_t> kv1(n_kmers), kv2(n_kmers);

  for(size_t i=0;i<nseqs;i++) {
    CStr seq1(intstr(s1[i].c_str()));
    len1 = s1[i].size();
    assign_kmer(kv1.data(), seq1.get(), kmer_size);
    CStr seq2(intstr(s2[i].c_str()));
    len2 = s2[i].size();
    assign_kmer(kv2.data(), seq2.get(), kmer_size);
    kdist[i] = ::kmer_dist(kv1.data(), len1, kv2.data(), len2, kmer_size);
  }
  return kdist;
}

std::vector<double> kord_dist(const std::vector<std::string> &s1, const std::vector<std::string> &s2,
                              int kmer_size, int SSE) {
  size_t len1 = 0, len2 = 0, maxlen=0;
  check_kmer_size(kmer_size, 1);

  size_t nseqs = s1.size();
  if(nseqs != s2.size()) { dada2_stop_invalid("Mismatched numbers of sequences."); }
  for(size_t i=0;i<nseqs;i++) {
    len1 = s1[i].size();
    len2 = s2[i].size();
    if(len1 > maxlen) { maxlen = len1; }
    if(len2 > maxlen) { maxlen = len2; }
  }

  std::vector<double> kdist(nseqs);
  std::vector<uint16_t> kord1(maxlen + 8), kord2(maxlen + 8);  // +8: SSE loads read whole blocks

  for(size_t i=0;i<nseqs;i++) {
    CStr seq1(intstr(s1[i].c_str()));
    len1 = s1[i].size();
    assign_kmer_order(kord1.data(), seq1.get(), kmer_size);
    CStr seq2(intstr(s2[i].c_str()));
    len2 = s2[i].size();
    assign_kmer_order(kord2.data(), seq2.get(), kmer_size);
    if(SSE==1) {
      kdist[i] = kord_dist_SSEi(kord1.data(), len1, kord2.data(), len2, kmer_size);
    } else {
      kdist[i] = ::kord_dist(kord1.data(), len1, kord2.data(), len2, kmer_size);
    }
  }
  return kdist;
}

std::vector<int> kmer_matches(const std::vector<std::string> &s1,
                              const std::vector<std::string> &s2, int kmer_size) {
  size_t len1 = 0, len2 = 0, maxlen=0;
  size_t klen1 = 0, klen2 = 0, klen_min = 0;
  check_kmer_size(kmer_size, 1);
  int matches=0;

  size_t nseqs = s1.size();
  if(nseqs != s2.size()) { dada2_stop_invalid("Mismatched numbers of sequences."); }
  for(size_t i=0;i<nseqs;i++) {
    len1 = s1[i].size();
    len2 = s2[i].size();
    if(len1 > maxlen) { maxlen = len1; }
    if(len2 > maxlen) { maxlen = len2; }
  }

  std::vector<int> kmatch(nseqs);
  std::vector<uint16_t> kord1(maxlen), kord2(maxlen);

  for(size_t i=0;i<nseqs;i++) {
    CStr seq1(intstr(s1[i].c_str()));
    len1 = s1[i].size();
    klen1 = len1 - kmer_size + 1;
    assign_kmer_order(kord1.data(), seq1.get(), kmer_size);
    CStr seq2(intstr(s2[i].c_str()));
    len2 = s2[i].size();
    klen2 = len2 - kmer_size + 1;
    assign_kmer_order(kord2.data(), seq2.get(), kmer_size);
    // Calculate matches
    matches=0;
    klen_min = klen1 < klen2 ? klen1 : klen2;
    for(size_t j=0;j<klen_min;j++) {
      if(kord1[j] == kord2[j]) { matches++; }
    }
    kmatch[i] = matches;
  }
  return kmatch;
}

std::vector<int> kdist_matches(const std::vector<std::string> &s1,
                               const std::vector<std::string> &s2, int kmer_size) {
  uint16_t dotsum = 0;
  check_kmer_size(kmer_size, 3);
  size_t n_kmers = (1 << (2*kmer_size));  // 4^k kmers

  size_t nseqs = s1.size();
  if(nseqs != s2.size()) { dada2_stop_invalid("Mismatched numbers of sequences."); }

  std::vector<int> kdist_match(nseqs);
  std::vector<uint16_t> kv1(n_kmers), kv2(n_kmers);

  for(size_t i=0;i<nseqs;i++) {
    CStr seq1(intstr(s1[i].c_str()));
    assign_kmer(kv1.data(), seq1.get(), kmer_size);
    CStr seq2(intstr(s2[i].c_str()));
    assign_kmer(kv2.data(), seq2.get(), kmer_size);
    // code from kmer_dist
    dotsum = 0;
    for(size_t j=0;j<n_kmers;j++) {
      dotsum += (kv1[j] < kv2[j] ? kv1[j] : kv2[j]);
    }
    kdist_match[i] = dotsum;
  }
  return kdist_match;
}

}  // namespace dada2
