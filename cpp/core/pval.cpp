// Adapted from DADA2 src/pval.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * Rcpp::ppois(n, lambda, lower.tail=false) -> ppois_upper (stats.h, Boost.Math based).
//  * Removed the commented-out copy of R's nmath ppois/pgamma (GPL) carried upstream.
//  * compute_lambda: single thread-safe implementation over an ErrMat view; NaN-safe range check.
//  * compute_lambda_packed: the same computation from cached (packed) substitutions.
/*
Pval.cpp contains the functions related to calculating the abundance pval in DADA2.
*/

#include "dada.h"
#include "stats.h"


/* b_p_update:
 Calculates the abundance p-value for each raw in the clustering.
 Depends on the lambda between the raw and its cluster, and the reads of each.
 */
void b_p_update(B *b, bool greedy, bool detect_singletons) {
  unsigned int i, r;
  Raw *raw;
  Bi *bi;
  double E_reads_center;
  for(i=0;i<b->nclust;i++) {
    bi = b->bi[i];
    if(bi->update_e) {
      for(r=0;r<bi->nraw;r++) {
        raw = bi->raw[r];
        raw->p = get_pA(raw, bi, detect_singletons);
      } // for(r=0;r<b->bi[i]->nraw;r++)
      bi->update_e = false;
    } // if(bi->update_e)
    
    if(greedy && bi->check_locks) { // Lock raw if its reads are less than the expected reads from just the center 
      for(r=0;r<bi->nraw;r++) {
        raw = bi->raw[r];
        E_reads_center = b->bi[i]->center->reads * raw->comp.lambda;
        if(E_reads_center > raw->reads) { raw->lock = true; }
///!          if(raw->lock) { dada2_log("Raw %i (%i reads): E_Center = %.2f, P = %.2f\n", raw->index, raw->reads, E_reads_center, raw->p); }
        if(raw == bi->center) { raw->lock = true; }
      }
      bi->check_locks = false; // Locking can only happen first time around
    } // if(greedy && bi->check_locks)
  } // for(i=0;i<b->nclust;i++)
}

// Calculate abundance pval for given reads and expected number of reads
// Pval is conditional on sequnce being present, unless prior evidence is true
double calc_pA(int reads, double E_reads, bool prior) {
  double norm, pval=1.;
  
  // Calculate pval from poisson cdf.
  // -1 since strict > being calculated, and want to include the observed count
  pval = ppois_upper(reads-1, E_reads);  // was Rcpp::ppois(n, E_reads, lower.tail=false): P(X > x)

  if(!prior) {
    // Calculate norm (since conditioning on sequence being present).
    norm = (1.0 - exp(-E_reads));
    if(norm < TAIL_APPROX_CUTOFF) {
      norm = E_reads - 0.5*E_reads*E_reads; 
      // Assumption: TAIL_APPROX_CUTOFF is small enough to terminate taylor expansion at 2nd order
    }
    pval = pval/norm;
  }
  
  return pval;
}

// Find abundance pval from a Raw in a Bi
double get_pA(Raw *raw, Bi *bi, bool detect_singletons) {
  unsigned int hamming;
  double lambda, E_reads, pval = 1.;
  
  lambda = raw->comp.lambda;
  hamming = raw->comp.hamming;

  if(raw->reads == 1 && !raw->prior && !detect_singletons) {   // Singleton. No abundance pval.
    pval=1.;
  } 
  else if(hamming == 0) { // Cluster center (or no mismatch to center)
    pval=1.;
  }
  else if(lambda == 0) { // Zero expected reads of this raw
    pval = 0.;
  }
  else { // Calculate abundance pval.
    // E_reads is the expected number of reads for this raw
    E_reads = lambda * bi->reads;
    pval = calc_pA(raw->reads, E_reads, raw->prior || detect_singletons);
  }
  return pval;
}

// This calculates lambda from a lookup table index by transition (row) and rounded quality (col)
// pydada2: the upstream serial (Rcpp-matrix) and thread-safe variants were identical apart from
// the matrix access; only the thread-safe one (with the quality-range check) is kept.
double compute_lambda(Raw *raw, Sub *sub, const ErrMat &errMat, bool use_quals) {
  const unsigned int ncol = errMat.ncol;
  const double *err_mat = errMat.data;
  int s, pos0, pos1, nti0, nti1, len1;
  double lambda;
  unsigned int tvec[SEQLEN];
  unsigned int qind[SEQLEN];
  
  if(!sub) { // NULL Sub, outside Kmer threshold
    return 0.0;
  }
  
  // Make vector that indexes as integers the transitions at each position in seq1
  // Index is 0: A->A, 1: A->C, ..., 4: C->A, ...
  len1 = raw->length;
  for(pos1=0;pos1<len1;pos1++) {
    nti1 = ((int) raw->seq[pos1]) - 1;
    if(nti1 == 0 || nti1 == 1 || nti1 == 2 || nti1 == 3) {
      tvec[pos1] = nti1*4 + nti1;
    } else {
      dada2_stop("Non-ACGT sequences in compute_lambda.");
    }
    if(use_quals) {
      // Turn quality into the index in the array
      qind[pos1] = raw->qual[pos1]; // qind = unsigned int
    } else {
      qind[pos1] = 0;
    }
    
    if( qind[pos1] > (ncol-1) ) {
      dada2_stop("Rounded quality exceeded range of err lookup table.");
    }
  }
  
  // Now fix the ones where subs occurred
  for(s=0;s<sub->nsubs;s++) {
    pos0 = sub->pos[s];
    if(pos0 < 0 || pos0 >= sub->len0) { dada2_stop("CL: Bad pos0: %i (len0=%i).", pos0, sub->len0); }
    pos1 = sub->map[sub->pos[s]];
    if(pos1 < 0 || pos1 >= len1) { dada2_stop("CL: Bad pos1: %i (len1=%i).", pos1, len1); }
    
    nti0 = ((int) sub->nt0[s]) - 1;
    nti1 = ((int) sub->nt1[s]) - 1;
    tvec[pos1] = nti0*4 + nti1;
  }
  
  // And calculate lambda
  lambda = 1.0;
  for(pos1=0;pos1<len1;pos1++) {
    lambda = lambda * err_mat[tvec[pos1]*ncol+qind[pos1]];
  }
  
  if(!(lambda >= 0 && lambda <= 1)) { dada2_stop("Bad lambda."); } // NaN-safe
  
  return lambda;
}

// pydada2: compute_lambda for substitutions packed as (pos1 << 16 | nti0*4+nti1) (comparison
// cache). Same checks, same transition vector and the same product order as compute_lambda.
double compute_lambda_packed(Raw *raw, const uint32_t *subs, unsigned int nsubs, const ErrMat &errMat, bool use_quals) {
  const unsigned int ncol = errMat.ncol;
  const double *err_mat = errMat.data;
  int pos1, nti1, len1;
  unsigned int s;
  double lambda;
  unsigned int tvec[SEQLEN];
  unsigned int qind[SEQLEN];

  len1 = raw->length;
  for(pos1=0;pos1<len1;pos1++) {
    nti1 = ((int) raw->seq[pos1]) - 1;
    if(nti1 == 0 || nti1 == 1 || nti1 == 2 || nti1 == 3) {
      tvec[pos1] = nti1*4 + nti1;
    } else {
      dada2_stop("Non-ACGT sequences in compute_lambda.");
    }
    if(use_quals) {
      qind[pos1] = raw->qual[pos1];
    } else {
      qind[pos1] = 0;
    }
    if( qind[pos1] > (ncol-1) ) {
      dada2_stop("Rounded quality exceeded range of err lookup table.");
    }
  }
  for(s=0;s<nsubs;s++) {
    pos1 = (int)(subs[s] >> 16);
    if(pos1 >= len1) { dada2_stop("CL: Bad pos1: %i (len1=%i).", pos1, len1); }
    tvec[pos1] = subs[s] & 0xFFFF;
  }
  lambda = 1.0;
  for(pos1=0;pos1<len1;pos1++) {
    lambda = lambda * err_mat[tvec[pos1]*ncol+qind[pos1]];
  }
  if(!(lambda >= 0 && lambda <= 1)) { dada2_stop("Bad lambda."); } // NaN-safe
  return lambda;
}
