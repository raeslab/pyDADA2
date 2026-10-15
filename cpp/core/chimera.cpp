// Adapted from DADA2 src/chimera.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * C_is_bimera -> dada2::is_bimera; BimeraTableParallel/C_table_bimera2 -> dada2::table_bimera
//    using dada2::parallel_for with a ROW-major (samples x sequences) table.
//  * Alignments are freed on every path (upstream leaked on the `continue` path of is_bimera).
//  * is_bimera_ptr: is_bimera on parents given by pointer (no copies of the parent strings).
#include <algorithm>

#include "api.h"
#include "dada.h"
#include "parallel.h"


namespace dada2 {

//------------------------------------------------------------------
// Determines whether sq is a perfect bimera of some combination from pars (was C_is_bimera).
//
// @param sq The query DNA sequence.
// @param pars Potential bimeric "parents".
// @param max_shift A \code{integer(1)} of the maximum alignment shift allowed.
bool is_bimera(const std::string &sq, const std::vector<std::string> &pars, bool allow_one_off,
               int min_one_off_par_dist, int match, int mismatch, int gap_p, int max_shift) {
  std::vector<const std::string *> ptrs(pars.size());
  for(size_t i=0;i<pars.size();i++) { ptrs[i] = &pars[i]; }
  return is_bimera_ptr(sq, ptrs, allow_one_off, min_one_off_par_dist, match, mismatch, gap_p, max_shift);
}

bool is_bimera_ptr(const std::string &sq, const std::vector<const std::string *> &pars, bool allow_one_off,
                   int min_one_off_par_dist, int match, int mismatch, int gap_p, int max_shift) {
  int left, right, left_oo, right_oo;
  char **al;
  int max_left=0, max_right=0;
  int oo_max_left=0, oo_max_right=0, oo_max_left_oo=0, oo_max_right_oo=0;
  bool rval = false;
  const int sqlen = sq.size();

  for(size_t i=0;i<pars.size() && rval==false;i++) {
    al = nwalign_vectorized2(sq.c_str(), sq.length(), pars[i]->c_str(), pars[i]->length(), (int16_t) match, (int16_t) mismatch, (int16_t) gap_p, 0, max_shift);  // Remember, alignments must be freed!
    get_lr(al, left, right, left_oo, right_oo, allow_one_off, max_shift);

    if((left+right) >= sqlen) { // Toss id/pure-shift/internal-indel "parents"
      al_free(al);
      continue;
    }
    if(left > max_left) { max_left=left; }
    if(right > max_right) { max_right=right; }

    // Need to evaluate whether parents are allowed for one-off models
    if(allow_one_off && get_ham_endsfree(al[0], al[1]) >= min_one_off_par_dist) {
      if(left > oo_max_left) { oo_max_left=left; }
      if(right > oo_max_right) { oo_max_right=right; }
      if(left_oo > oo_max_left_oo) { oo_max_left_oo=left_oo; }
      if(right_oo > oo_max_right_oo) { oo_max_right_oo=right_oo; }
    }

    // Found bimeric model?
    if((max_right+max_left)>=sqlen) {
      rval = true;
    }
    if(allow_one_off) {
      if((oo_max_left+oo_max_right_oo)>=sqlen || (oo_max_left_oo+oo_max_right)>=sqlen) {
        rval=true;
      }
    }
    al_free(al);
  }

  return(rval);
}

// Flags bimeras in each sample of a sequence table (was BimeraTableParallel + C_table_bimera2).
// mat: nrow (samples) x ncol (sequences), ROW-major (upstream: R column-major matrix).
std::pair<std::vector<int>, std::vector<int>> table_bimera(
    const int32_t *mat, int nrow, int ncol, const std::vector<std::string> &seqs, double min_fold,
    int min_abund, bool allow_one_off, int min_one_off_par_dist, int match, int mismatch,
    int gap_p, int max_shift, int n_threads) {
  if((int)seqs.size() != ncol) { dada2_stop_invalid("Number of sequences does not match the number of table columns."); }
  std::vector<int> C_flags(ncol, 0);
  std::vector<int> C_sams(ncol, 0);
  // mat(i,j) --> vals[i*ncol+j]
  const int32_t *vals = mat;

  auto body = [&](std::size_t begin, std::size_t end) {
    int i,k,nsam,nflag,sqlen,left,right,left_oo,right_oo,max_left,max_right;
    int oo_max_left, oo_max_right, oo_max_left_oo, oo_max_right_oo;
    char **al;
    std::vector<int> lefts(ncol);
    std::vector<int> rights(ncol);
    std::vector<int> lefts_oo(ncol);
    std::vector<int> rights_oo(ncol);
    std::vector<bool> allowed(ncol);

    for(std::size_t j=begin;j<end;j++) {
      nsam=0; nflag=0;
      sqlen = seqs[j].size();
      std::fill(lefts.begin(), lefts.end(), -1);
      std::fill(rights.begin(),rights.end(),-1);
      if(allow_one_off) {
        std::fill(lefts_oo.begin(), lefts_oo.end(), -1);
        std::fill(rights_oo.begin(),rights_oo.end(),-1);
        std::fill(allowed.begin(), allowed.end(), false);
      }
      for(i=0;i<nrow;i++) { // Evaluate in each sample (row)
        if(vals[(size_t)i*ncol+j]<=0) { continue; }
        nsam++;
        max_left=0; max_right=0;
        oo_max_left=0; oo_max_right=0; oo_max_left_oo=0; oo_max_right_oo=0;
        for(k=0;k<ncol;k++) { // Compare with all possible parents
          if(vals[(size_t)i*ncol+k]>(min_fold*vals[(size_t)i*ncol+j]) && vals[(size_t)i*ncol+k]>=min_abund) {
            if(lefts[k]<0) { // Comparison not yet done to this potential parent
              al = nwalign_vectorized2(seqs[j].c_str(), seqs[j].length(), seqs[k].c_str(), seqs[k].length(), (int16_t) match, (int16_t) mismatch, (int16_t) gap_p, 0, max_shift);  // Remember, alignments must be freed!
              get_lr(al, left, right, left_oo, right_oo, allow_one_off, max_shift);
              if(allow_one_off && get_ham_endsfree(al[0], al[1]) >= min_one_off_par_dist) {
                allowed[k]=true;
              }

              if((left+right) < sqlen) {
                lefts[k]=left;
                rights[k]=right;
                if(allow_one_off) {
                  lefts_oo[k] = left_oo;
                  rights_oo[k] = right_oo;
                }
              } else {  // Ignore id/pure-shift/internal-indel "parents"
                lefts[k]=0;
                rights[k]=0;
                if(allow_one_off) {
                  lefts_oo[k] = 0;
                  rights_oo[k] = 0;
                }
              }
              al_free(al);
            }
            // Now compare to best parents yet found
            if(lefts[k] > max_left) { max_left=lefts[k]; }
            if(rights[k] > max_right) { max_right=rights[k]; }
            if(allow_one_off && allowed[k]) {
              if(lefts[k] > oo_max_left) { oo_max_left=lefts[k]; }
              if(rights[k] > oo_max_right) { oo_max_right=rights[k]; }
              if(lefts_oo[k] > oo_max_left_oo) { oo_max_left_oo=lefts_oo[k]; }
              if(rights_oo[k] > oo_max_right_oo) { oo_max_right_oo=rights_oo[k]; }
            }
          } // if(vals[i+k*nrow]>vals[i+j*nrow] && vals[i+k*nrow]>=2)
        } // for(k=0;k<mat.ncol();k++)

        // Flag if chimeric model exists
        if((max_right+max_left)>=sqlen) {
          nflag++;
        } else if(allow_one_off) {
          if((oo_max_left+oo_max_right_oo)>=sqlen || (oo_max_left_oo+oo_max_right)>=sqlen) {
            nflag++;
          }
        }
      } // for(i=0;i<mat.nrow();i++)

      C_flags[j] = nflag;
      C_sams[j] = nsam;
    } // for(std::size_t j=begin;j<end;j++)
  };
  dada2::parallel_for(0, ncol, body, 1, n_threads);

  return {C_flags, C_sams};
}

}  // namespace dada2

// Internal function to get hamming distance between aligned seqs
//  without counting end gaps
int get_ham_endsfree(const char *seq1, const char *seq2) {
  int i, j, pos, ham;
  size_t len=strlen(seq2);
  bool gap1, gap2;
  // Find start of internal part of alignment
  i=0;
  gap1 = (seq1[i]=='-'); 
  gap2 = (seq2[i]=='-');
  while(gap1 || gap2) {
    i++;
    gap1 = (gap1 && (seq1[i]=='-'));
    gap2 = (gap2 && (seq2[i]=='-'));
  }
  // Find end of internal part of alignment
  j=len-1;
  gap1 = (seq1[j]=='-'); 
  gap2 = (seq2[j]=='-');
  while(gap1 || gap2) {
    j--;
    gap1 = (gap1 && (seq1[j]=='-'));
    gap2 = (gap2 && (seq2[j]=='-'));
  }
  // Calculate hamming distance over internal part
  ham=0;
  for(pos=i;pos<=j;pos++) {
    if(seq1[pos] != seq2[pos]) { ham++; }
  }
  return(ham);
}

// Internal function to get left and right overlaps, 
//   with or without allowing one mismatch
void get_lr(char **al, int &left, int &right, int &left_oo, int &right_oo, bool allow_one_off, int max_shift) {
  size_t len = strlen(al[0]);
  int pos=0; left=0;
  while(al[0][pos] == '-' && pos<len) {
    pos++; // Scan in until query starts
  }
  while(al[1][pos] == '-' && pos<max_shift) {
    pos++; left++; // Credit as ends-free coverage until parent starts
  }
  while(pos<len && al[0][pos] == al[1][pos]) {
    pos++; left++; // Credit as covered until a mismatch
  }
  if(allow_one_off) {
    // Step forward, and credit to one-off further matches (and this mismatch if not a gap)
    left_oo = left;
    pos++;
    if(pos<len && al[0][pos] != '-') { left_oo++; }
    while(pos<len && al[0][pos] == al[1][pos]) {
      pos++; left_oo++;
    }
  }
  
  pos=len-1; right=0;
  while(al[0][pos] == '-' && pos >= 0) { 
    pos--;
  }
  while(al[1][pos] == '-' && pos>+(len-max_shift)) {
    pos--; right++;
  }
  while(pos>=0 && al[0][pos] == al[1][pos]) {
    pos--; right++;
  }
  if(allow_one_off) {
    // Step forward, and credit to one-off further matches (and this mismatch if not a gap)
    right_oo = right;
    pos--;
    if(pos>=0 && al[0][pos] != '-') { right_oo++; }
    while(pos>=0 && al[0][pos] == al[1][pos]) {
      pos--; right_oo++;
    }
  }
}
