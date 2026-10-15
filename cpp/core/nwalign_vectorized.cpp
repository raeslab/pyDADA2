// Adapted from DADA2 src/nwalign_vectorized.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * C_nwvec -> dada2::nwvec (api.h); the R band argument was cast through (size_t), now passed as int (same value).
//  * __restrict__ -> DADA2_RESTRICT (MSVC); debug printer removed; Rcpp::stop -> dada2_stop.
//  * nwalign_vectorized2: the diagonal scores are filled from a reversed copy of s1, so both
//    sequences are read contiguously and the loop vectorises (same values); the DP matrices and
//    traceback buffers are thread-local and reused between calls instead of malloc'ed per call
//    (every cell that is read is written first in each call, checked with valgrind memcheck).

#include <vector>

#include "api.h"
#include "dada.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

// Precedence is up>left>diag
void dploop_vec(int16_t *DADA2_RESTRICT ptr_left, int16_t *DADA2_RESTRICT ptr_diag, int16_t *DADA2_RESTRICT ptr_up, int16_t *DADA2_RESTRICT d, int16_t *DADA2_RESTRICT p, int16_t gap_p, size_t n) {
  int16_t left, diag, up, entry, pentry;
  size_t i = 0;
  
  while(i<n) {
    left = *ptr_left + gap_p;
    diag = *ptr_diag;
    up = *ptr_up + gap_p;
    
    entry = up >= left ? up : left;
    pentry = up >= left ? 3 : 2;
    pentry = entry >= diag ? pentry : 1;
    entry = entry >= diag ? entry : diag;
    
    *d = entry; // Vectorizes if this points to another array
    *p = pentry;
    
    ptr_left++;
    ptr_diag++;
    ptr_up++;
    d++;
    p++;
    i++;
  }
}

// Precedence is left>up>diag
void dploop_vec_swap(int16_t *DADA2_RESTRICT ptr_left, int16_t *DADA2_RESTRICT ptr_diag, int16_t *DADA2_RESTRICT ptr_up, int16_t *DADA2_RESTRICT d, int16_t *DADA2_RESTRICT p, int16_t gap_p, size_t n) {
  int16_t left, diag, up, entry, pentry;
  size_t i = 0;
  
  while(i<n) {
    left = *ptr_left + gap_p;
    diag = *ptr_diag;
    up = *ptr_up + gap_p;
    
    entry = left >= up ? left : up;
    pentry = left >= up ? 2 : 3;
    pentry = entry >= diag ? pentry : 1;
    entry = entry >= diag ? entry : diag;
    
    *d = entry; // Vectorizes if this points to another array
    *p = pentry;
    
    ptr_left++;
    ptr_diag++;
    ptr_up++;
    d++;
    p++;
    i++;
  }
}

namespace {

// Per-thread scratch buffers of nwalign_vectorized2 (grown as needed, released when very large).
struct NwScratch {
  std::vector<int16_t> d, p, diag_buf;
  std::vector<char> s1r, al0, al1;
  void trim() {
    const size_t big = (size_t)16 << 20;  // keep up to 16 MiB per buffer between calls
    if(d.capacity() * sizeof(int16_t) > big) { std::vector<int16_t>().swap(d); std::vector<int16_t>().swap(p); }
  }
};

NwScratch &nw_scratch() {
  thread_local NwScratch scratch;
  return scratch;
}

}  // namespace

char **nwalign_vectorized2(const char *s1, size_t len1, const char *s2, size_t len2, int16_t match, int16_t mismatch, int16_t gap_p, int16_t end_gap_p, int band) {
  size_t row, col, ncol, nrow, foo;
  size_t i,j;
  int16_t d_free;
  size_t start_col;
  size_t col_min, col_max, even;
  size_t i_max, j_min;
  int16_t *ptr_left, *ptr_diag, *ptr_up, *ptr_d, *ptr_p;
  bool swap = false;
  bool recalc_left = false, recalc_right = false;
  const char *ptr_const_char;
  char *ptr_char;

  if(len1 > len2) { // Ensure s1 is the shorter sequences
    ptr_const_char = s1;
    s1 = s2;
    s2 = ptr_const_char;
    swap = true;
    foo = len1;
    len1 = len2;
    len2 = foo;
  }
  if(band < 0) { band = len2; }
  
  // Allocate the DP matrices
  start_col = 1 + (1+(band<len1 ? band : len1))/2;
  //  ncol = 3 + (len1+len2+1)/2; // 3 = left boundary + center + right boundary !!!
  ncol = 2 + start_col + ((len2-len1+band)<len2 ? (len2-len1+band) : len2)/2;
  nrow = len1 + len2 + 1;
  NwScratch &scratch = nw_scratch();
  if(scratch.d.size() < ncol * nrow) { scratch.d.resize(ncol * nrow); scratch.p.resize(ncol * nrow); }
  if(scratch.diag_buf.size() < ncol) scratch.diag_buf.resize(ncol);
  if(scratch.s1r.size() < len1) scratch.s1r.resize(len1);
  if(scratch.al0.size() < nrow + 1) { scratch.al0.resize(nrow + 1); scratch.al1.resize(nrow + 1); }
  int16_t *d = scratch.d.data();
  int16_t *p = scratch.p.data();
  int16_t *diag_buf = scratch.diag_buf.data();
  char *s1r = scratch.s1r.data();  // s1 reversed: s1[i] == s1r[len1-1-i]
  for(i=0;i<len1;i++) s1r[i] = s1[len1-1-i];
  
  // For banding issues later on
  int16_t fill_val = INT16_MIN - MIN(MIN(mismatch, gap_p), MIN(match, 0));
  for(row=0;row<nrow;row++) {
    d[row*ncol] = fill_val;
    d[row*ncol+1] = fill_val;
    d[row*ncol+ncol-2] = fill_val;
    d[row*ncol+ncol-1] = fill_val;
  }

  // Fill out starting point
  d[start_col] = 0;
  p[start_col] = 0; // Should never be queried
  
  // Fill out "left" "column" of d, p.
  row=1;
  col=start_col-1;
  ptr_d=&d[ncol]; // start of 1st row
  ptr_p=&p[ncol];
  d_free = end_gap_p; // end-gap score
  while(row < (1 + (band < len1 ? band : len1))) {
    ptr_d[col] = d_free; // end gap
    ptr_p[col] = 3;
    if(row%2==0) {
      col--;
    }
    row++;
    d_free += end_gap_p;
    ptr_d += ncol;
    ptr_p += ncol;
  }
  
  // Fill out "top" "row" of d, p.
  row=1;
  col=start_col;
  ptr_d=&d[ncol]; // start of 1st row
  ptr_p=&p[ncol];
  d_free = end_gap_p; // end-gap score
  while(row < (1+(band+len2-len1 < len2 ? band+len2-len1 : len2))) {
    ptr_d[col] = d_free;
    ptr_p[col] = 2;
    if(row%2 == 1) {
      col++;
    }
    row++;
    d_free += end_gap_p;
    ptr_d += ncol;
    ptr_p += ncol;
  }

  // Fill out DP matrix (Row 0/1 taken care of by ends-free)
  row = 2;
  col_min = start_col; // Do not fill out the ends-free cells
  col_max = start_col;
  i_max = 0; // 1st nt
  j_min = 0; // 1st nt
  even = 1; // TRUE
  
  while(row <= (len1+len2)) {
      // Fill out row
    // == for(col=col_min,i=i_max,j=j_min;col<=col_max;col++,i--,j++)
    //      diag_buf[col] = d[(row-2)*ncol + col] + (s1[i] == s2[j] ? match : mismatch);
    {
      const char *DADA2_RESTRICT a = s1r + (len1 - 1 - i_max);  // a[t] == s1[i_max - t]
      const char *DADA2_RESTRICT b = s2 + j_min;
      const int16_t *DADA2_RESTRICT dprev = &d[(row-2)*ncol + col_min];
      int16_t *DADA2_RESTRICT dbuf = &diag_buf[col_min];
      const size_t nn = col_max + 1 - col_min;
      for(size_t t=0;t<nn;t++) dbuf[t] = dprev[t] + (a[t] == b[t] ? match : mismatch);
    }
    ptr_left = &d[(row-1)*ncol + col_min-even];
    ptr_diag = &diag_buf[col_min];
    ptr_up = &d[(row-1)*ncol + col_min+1-even];
    if(swap) {
      dploop_vec_swap(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    } else {
      dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    }
    
    // Offset to the boundary after top wedge (w/ prefilled boundary) is done
    if(row==(band<len1 ? band : len1)) { 
      col_min--;
      i_max++;
      j_min--;
    }
    if(row == (band+len2-len1 < len2 ? band+len2-len1 : len2)) {
      col_max++;
    }
    
    // Recalculate ends-free cells for lower boundary
    if(end_gap_p > gap_p) {
      // Left column
      if(recalc_left) { // past first row of lower tri
        d_free = ptr_left[0] + end_gap_p; // first cell of left array
        if(d_free > d[row*ncol + col_min]) { // ends-free gap is better
          d[row*ncol + col_min] = d_free;
          p[row*ncol + col_min] = 2;
        } else if(!swap && d_free == d[row*ncol + col_min] && p[row*ncol + col_min] == 1) { // left gap takes precedence over diagonal move (for consistency)
          p[row*ncol + col_min] = 2;
        } else if(swap && d_free == d[row*ncol + col_min] && p[row*ncol + col_min] != 2) { // left-is-up, and takes precedence other moves
          p[row*ncol + col_min] = 2;
        }
      }
      if(i_max == len1-1) { recalc_left = true; }
      // Right column
      if(recalc_right) {
        d_free = ptr_up[col_max-col_min] + end_gap_p; // last cell of up array
        if(d_free > d[row*ncol + col_max]) { // ends-free gap is better
          d[row*ncol + col_max] = d_free;
          p[row*ncol + col_max] = 3;
        } else if(!swap && d_free == d[row*ncol + col_max] && p[row*ncol + col_max] != 3) { // up gap takes precedence over left or diagonal move (for consistency)
          p[row*ncol + col_max] = 3;
        } else if(swap && d_free == d[row*ncol + col_max] && p[row*ncol + col_max] == 1) { // up-is-left, and takes precedence over diagonal
          p[row*ncol + col_max] = 3;
        }
      }
      if((row+1)/2 + col_max - start_col == len2) { recalc_right = true; }
      // j_max (1-index) = (row+1)/2 + col_max - start_col
    }
    // Update the indices
    if(row < band && row < len1) { // upper tri for seq1
      if(even) { col_min--; }
      i_max++;
    } else if(i_max < len1-1) { // banded area
      if(band%2 == 0) {
        if(even) { j_min++; }
        else { i_max++; }
      } else { // odd band
        if(even) { col_min--; i_max++; }
        else { col_min++; j_min++; }
      }
    } else { // lower tri for seq1
      if(!even) { col_min++; }
      j_min++;
    }

    if(row<(band+len2-len1 < len2 ? band+len2-len1 : len2)) {
      if(!even) { col_max++; }
    } else if((row+1)/2 + col_max - start_col < len2) { // "j_max" (1-index) < len2
      if((band+len2-len1) % 2 == 0) { // even band (including the extra band from length difference)
        if(even) { col_max--; }
        else { col_max++; }
      } // no action on odd band
    } else {
      if(even) { col_max--; }
    }
    
    row++;
    even = 1 - even;
  }

//  if(debug) {
//    dada2_log("D MATRIX:.....\n");
//    parr(d, len1+len2+1, ncol);
//    dada2_log("\n");
//    dada2_log("Score: %d\n", d[(len1+len2)*ncol + (2*start_col+len2-len1)/2]);
//  }

  char *al0 = scratch.al0.data();
  char *al1 = scratch.al1.data();
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;
  
  while ( i > 0 || j > 0 ) {
    switch ( p[(i+j)*ncol + (2*start_col+j-i)/2] ) {
      case 1:
        al0[len_al] = s1[--i];
        al1[len_al] = s2[--j];
        break;
      case 2:
        al0[len_al] = '-';
        al1[len_al] = s2[--j];
        break;
      case 3:
        al0[len_al] = s1[--i];
        al1[len_al] = '-';
        break;
      default:
        dada2_stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  
  // Return to input ordering
  if(swap) {
    ptr_char = al0;
    al0 = al1;
    al1 = ptr_char;
  }
  
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  dada2_stop("Failed memory allocation.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  dada2_stop("Failed memory allocation.");
  
  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  scratch.trim();
  return al;
}

namespace dada2 {

// Vectorized Needleman-Wunsch alignment of s1[i] with s2[i] (was C_nwvec).
std::vector<std::pair<std::string, std::string>> nwvec(const std::vector<std::string> &s1,
                                                       const std::vector<std::string> &s2,
                                                       int match, int mismatch, int gap_p,
                                                       int band, bool endsfree) {
  char **al;
  if(s1.size() != s2.size()) {
    dada2_stop_invalid("Character vectors to be aligned must be of equal length.");
  }
  std::vector<std::pair<std::string, std::string>> rval(s1.size());

  for(size_t i=0;i<s1.size();i++) {
    if(endsfree) {
      al = nwalign_vectorized2(s1[i].c_str(), s1[i].length(), s2[i].c_str(), s2[i].length(), (int16_t) match, (int16_t) mismatch, (int16_t) gap_p, 0, band);
    } else {
      al = nwalign_vectorized2(s1[i].c_str(), s1[i].length(), s2[i].c_str(), s2[i].length(), (int16_t) match, (int16_t) mismatch, (int16_t) gap_p, (int16_t) gap_p, band);
    }
    rval[i] = {std::string(al[0]), std::string(al[1])};
    al_free(al);
  }
  return rval;
}

}  // namespace dada2
