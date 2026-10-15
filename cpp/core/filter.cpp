// Adapted from DADA2 src/filter.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * C_matchRef -> dada2::match_ref, C_matrixEE -> dada2::matrix_ee (row-major input, NA_QUAL for NA).

#include "api.h"
#include "dada.h"

namespace dada2 {

// Counts word_size-mers of each sequence found in ref (circular), was C_matchRef.
std::vector<int> match_ref(const std::vector<std::string> &seqs, const std::string &ref_in,
                           int word_size_in, bool non_overlapping) {
  int i,j;
  if(word_size_in < 1) { dada2_stop_invalid("word_size must be positive."); }
  unsigned int word_size = word_size_in;
  std::string ref = ref_in;
  std::unordered_set<std::string> phash; ///!
  std::vector<int> rval(seqs.size(), 0);

  unsigned int len = ref.size();
  ref.append(ref, 0, word_size);

  for(i=0;i<(int)len;i++) {
    phash.insert(ref.substr(i,word_size));
  }

  for(i=0;i<(int)seqs.size();i++) {
    len=seqs[i].size();
    if(len<word_size) { continue; }

    for(j=0;j<=(int)(len-word_size);j++) {
      if(phash.count(seqs[i].substr(j, word_size))) {
        rval[i]++;
        if(non_overlapping) { j+=word_size; }
      }
    }
  }
  return(rval);
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
