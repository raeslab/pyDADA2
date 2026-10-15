// Adapted from DADA2 src/cluster.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * b_compare/b_compare_parallel merged into one b_compare using dada2::parallel_for.
//  * Rcpp::NumericMatrix -> ErrMat; Rcpp::stop/Rprintf -> dada2_stop/dada2_log.
//  * b_shuffle2 scratch arrays are std::vectors (exception safety).
//  * b_bud sets birth_from for prior-born ("P") clusters (upstream left it unset: garbage).
//  * Lambda range checks are NaN-safe.
//  * b_compare can use a cache of earlier comparisons with the same centre (see CentreCache in
//    dada.h): the k-mer screen decision is re-applied to the cached k-mer distance and lambda is
//    recomputed from the cached substitutions exactly as compute_lambda does.
#include "dada.h"
#include "parallel.h"

/********* ALGORITHM LOGIC *********/

/*
 compare:
Performs alignments and computes lambda for all raws to the specified Bi
Stores only those that can possibly be recruited to this Bi

pydada2: upstream had a serial b_compare and an RcppParallel b_compare_parallel that produced
identical comparisons. They are merged: comparisons are computed with dada2::parallel_for
(serially when n_threads <= 1) into a per-raw buffer, then stored sequentially in raw order
exactly as upstream's parallel variant did.
*/
void b_compare(B *b, unsigned int i, const ErrMat &errMat,
               int match, int mismatch, int gap_pen, int homo_gap_pen,
               bool use_kmers, double kdist_cutoff, int band_size, bool vectorized_alignment,
               int SSE, bool gapless, bool greedy, bool verbose, int n_threads) {
  unsigned int index;
  double lambda;
  Raw *raw;
  Comparison comp;

  if(verbose) { dada2_log("C%iLU:", i); }

  // Perform all comparisons (thread-safe: each index writes only its own slot)
  std::vector<Comparison> comps(b->nraw);
  std::vector<unsigned char> aligned(b->nraw, 0), shrouded(b->nraw, 0);
  Raw *center = b->bi[i]->center;
  const bool use_quals = b->use_quals;

  // pydada2: comparison cache (alignments and k-mer distances do not depend on the error matrix).
  // New results are collected per raw and merged into the cache serially below.
  CentreCache *cc = b->cache ? cache_centre(b->cache, center->index, b->nraw) : NULL;
  struct NewSubs {  // owns the subs kept for the cache until the merge (exception safety)
    std::vector<Sub *> v;
    ~NewSubs() { for(Sub *s : v) sub_free(s); }
  } new_subs;
  std::vector<double> new_kdist;
  std::vector<uint8_t> new_state;  // 0: nothing new, 1: screened, 2: aligned
  if(cc) {
    new_subs.v.assign(b->nraw, NULL);
    new_kdist.assign(b->nraw, 0.0);
    new_state.assign(b->nraw, 0);
  }

  auto body = [&](std::size_t begin, std::size_t end) {
    Raw *raw;
    Sub *sub;
    for(std::size_t index=begin;index<end;index++) {
      raw = b->raw[index];
      comps[index].i = i;
      comps[index].index = index;
      if(greedy && (raw->reads > center->reads)) {
        sub = NULL;
      } else if(greedy && raw->lock) { // was (raw->p > 0.5)
        sub = NULL;
      } else if(cc && cc->nsubs[index] != CACHE_UNKNOWN &&
                (cc->nsubs[index] != CACHE_SCREENED || (use_kmers && cc->kdist[index] > kdist_cutoff))) {
        // Cached: same outcome as sub_new + compute_lambda below
        aligned[index] = 1;
        if(use_kmers && cc->kdist[index] > kdist_cutoff) { // k-mer screen (raw_align)
          shrouded[index] = 1;
          comps[index].lambda = 0.0; // compute_lambda(raw, NULL, ...)
          comps[index].hamming = -1;
        } else {
          const uint32_t *packed = cc->subs.data() + cc->off[index];
          comps[index].lambda = b->lam ? compute_lambda_tab_packed(*b->lam, raw, packed, cc->nsubs[index])
                                       : compute_lambda_packed(raw, packed, cc->nsubs[index], errMat, use_quals);
          comps[index].hamming = cc->nsubs[index];
        }
        continue;
      } else {
        double kdist = 0.0;
        sub = sub_new(center, raw, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless, cc ? &kdist : NULL);
        aligned[index] = 1;
        if(!sub) { shrouded[index] = 1; }
        if(cc) {
          new_kdist[index] = kdist;
          new_state[index] = sub ? 2 : 1;
        }
      }

      // Make comparison object
      try {
        comps[index].lambda = b->lam ? compute_lambda_tab(*b->lam, raw, sub)
                                     : compute_lambda(raw, sub, errMat, use_quals);
      } catch(...) {
        sub_free(sub);
        throw;
      }
      if(sub) {
        comps[index].hamming = sub->nsubs;
      } else {
        comps[index].hamming = -1;
      }
      if(cc && sub) { new_subs.v[index] = sub; } // kept for the cache
      else { sub_free(sub); }
    }
  };
  dada2::parallel_for(0, b->nraw, body, GRAIN_SIZE, n_threads);

  if(cc) { // store the new comparisons (serial, raw order)
    for(index=0; index<b->nraw; index++) {
      if(new_state[index] == 1) {
        cc->kdist[index] = new_kdist[index];
        cc->nsubs[index] = CACHE_SCREENED;
      } else if(new_state[index] == 2) {
        Sub *sub = new_subs.v[index];
        if(cache_reserve(b->cache, sub->nsubs * sizeof(uint32_t))) {
          cc->kdist[index] = new_kdist[index];
          cc->off[index] = (uint32_t)cc->subs.size();
          for(unsigned int s=0;s<sub->nsubs;s++) {
            uint32_t pos1 = sub->map[sub->pos[s]];
            uint32_t t = (uint32_t)((sub->nt0[s]-1)*4 + (sub->nt1[s]-1));
            cc->subs.push_back((pos1 << 16) | t);
          }
          cc->nsubs[index] = (uint16_t)sub->nsubs;
        }
        sub_free(sub);
        new_subs.v[index] = NULL;
      }
    }
  }

  // Selectively store
  for(index=0; index<b->nraw; index++) {
    b->nalign += aligned[index];
    b->nshroud += shrouded[index];
    raw = b->raw[index];
    comp = comps[index];
    lambda = comp.lambda;
    if(!(lambda >= 0 && lambda <= 1)) dada2_stop("Lambda out-of-range error."); // NaN-safe

    // Store self-lambda
    if(index == b->bi[i]->center->index) {
      b->bi[i]->self = lambda;
    }

    // Store comparison if potentially useful
    if(lambda * b->reads > raw->E_minmax) { // This cluster could attract this raw
      if(lambda * b->bi[i]->center->reads > raw->E_minmax) { // Better E_minmax, set
        raw->E_minmax = lambda * b->bi[i]->center->reads;
      }
      b->bi[i]->comp.push_back(comp);
      if(i==0 || raw == b->bi[i]->center) { // Update on init (i=0) or if the center (as b_bud doesn't update raw->comp)
        raw->comp = comp;
      } // Handle init better?
    }
  }
}

/* b_shuffle2:
 move each sequence to the bi that produces the highest expected
 number of that sequence. The center of a Bi cannot leave.
*/
bool b_shuffle2(B *b) {
  unsigned int i, cind, index;
  double e;
  bool shuffled = false;
  Comparison *comp;
  Raw *raw;
  
  std::vector<double> emax(b->nraw);
  std::vector<Comparison *> compmax(b->nraw);

  // Initialize emax/imax off of cluster 0
  // Comparisons to all raws exist in cluster 0, in index order
  for(index=0;index<b->nraw;index++) {
    compmax[index] = &b->bi[0]->comp[index];
    emax[index] = compmax[index]->lambda * b->bi[0]->reads;
  }
  
  // Iterate over comparisons, find comparison with best E for each raw
  for(i=1;i<b->nclust;i++) {
    for(cind=0;cind<b->bi[i]->comp.size();cind++) {
      comp = &b->bi[i]->comp[cind];
      index = comp->index;
      e = comp->lambda * b->bi[i]->reads;
      if(e > emax[index]) { // better E
        compmax[index] = comp;
        emax[index] = e;
      }
    }
  }
  
  // Iterate over raws, if best i different than current, move
  for(i=0;i<b->nclust;i++) {
    // IMPORTANT TO ITERATE BACKWARDS DUE TO BI_POP_RAW!!!!!!
    for(int r=b->bi[i]->nraw-1; r>=0; r--) {
      raw = b->bi[i]->raw[r];
      // If a better cluster was found, move the raw to the new bi
      if(compmax[raw->index]->i != i) {
        if(raw->index == b->bi[i]->center->index) {  // Check if center
          if(VERBOSE) { dada2_log("Warning: Shuffle blocked the center of a Bi from leaving."); }
          continue;
        }
        // Moving raw
        bi_pop_raw(b->bi[i], r);
        bi_add_raw(b->bi[compmax[raw->index]->i], raw);
        // Assign raw the Comparison of its new cluster
        raw->comp = *compmax[raw->index];
        shuffled = true;  
      }  
    } // for(r=0;r<b->bi[i]->nraw;r++)
  }

  return shuffled;
}

/* b_bud:
 Finds the minimum p-value. If significant, creates a new cluster and moves the
 raws from the raw with the minimum p-value to the new cluster.
 Returns index of new cluster, or 0 if no new cluster added.
*/

int b_bud(B *b, double min_fold, int min_hamming, int min_abund, bool verbose) {
  int i, r, hamming;
  int mini = -1, minr = -1, mini_prior = -1, minr_prior = -1; // Negative initialization
  double pA, pP;
  double lambda, expected;
  Raw *raw, *minraw, *minraw_prior;
  minraw = b->bi[0]->center; // Assumes that complete alignment/pval calcs were performed in the init cluster
  minraw_prior = b->bi[0]->center; // Assumes that complete alignment/pval calcs were performed in the init cluster
  
  // Find i, r indices and value of minimum pval.
  for(i=0;i<b->nclust;i++) {
    for(r=1; r<b->bi[i]->nraw; r++) { // r=0 is the center
      raw = b->bi[i]->raw[r];
///      if(b->bi[i]->center->index == raw->index) { continue; } // Don't bud centers
      if(raw->reads < min_abund) { continue; }
      hamming = raw->comp.hamming;
      lambda = raw->comp.lambda;

      // Calculate the fold over-abundance and the hamming distance to this raw
      if(hamming >= min_hamming) { // Only those passing the hamming/fold screens can be budded
        if(min_fold <= 1 || ((double) raw->reads) >= min_fold * lambda * b->bi[i]->reads) {  
          if((raw->p < minraw->p) ||
             ((raw->p == minraw->p && raw->reads > minraw->reads))) { // Most significant
            mini = i; minr = r;
            minraw = raw;
          }
          if(raw->prior && ((raw->p < minraw_prior->p) ||
             (raw->p == minraw_prior->p && raw->reads > minraw_prior->reads))) { // Most significant
            mini_prior = i; minr_prior = r;
            minraw_prior = raw;
          }
        }
      }
    }
  }

///!  dada2_log("Prior... index: %i, i:%i, r:%i, p:%.2e\n", minraw_prior->index, mini_prior, minr_prior, minraw_prior->p); ///!
  // Bonferoni correct the abundance pval by the number of raws and compare to OmegaA
  // (quite conservative, although probably unimportant given the abundance model issues)
  pA = minraw->p * b->nraw;
  pP = minraw_prior->p;
  if(pA < b->omegaA && mini >= 0) {  // A significant abundance pval
    expected = minraw->comp.lambda * b->bi[mini]->reads;
    raw = bi_pop_raw(b->bi[mini], minr);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "A");
    b->bi[i]->birth_from = mini;
    b->bi[i]->birth_pval = pA;
    b->bi[i]->birth_fold = raw->reads/expected;
    b->bi[i]->birth_e = expected;
    b->bi[i]->birth_comp = minraw->comp;
    
    // Add raw to new cluster.
    bi_add_raw(b->bi[i], raw);
    bi_assign_center(b->bi[i]);
    if(verbose) { dada2_log(", Division (naive): Raw %i from Bi %i, pA=%.2e", raw->index, mini, pA); }
    return i;
  } else if (pP < b->omegaP && mini_prior >= 0) {  // A significant prior-abundance pval
    expected = minraw_prior->comp.lambda * b->bi[mini_prior]->reads;
    raw = bi_pop_raw(b->bi[mini_prior], minr_prior);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "P");
    b->bi[i]->birth_from = mini_prior; // pydada2: upstream left birth_from unset (garbage) here
    b->bi[i]->birth_pval = pP;
    b->bi[i]->birth_fold = raw->reads/expected;
    b->bi[i]->birth_e = expected;
    b->bi[i]->birth_comp = minraw_prior->comp;
    
    // Add raw to new cluster.
    bi_add_raw(b->bi[i], raw);
    bi_assign_center(b->bi[i]);
    if(verbose) { dada2_log(", Division (prior): Raw %i from Bi %i, pP=%.2e", raw->index, mini_prior, pP); }
    return i;
  }

  if(verbose) { dada2_log(", No Division. Minimum pA=%.2e (Raw %i w/ %i reads in Bi %i).", pA, minraw->index, minraw->reads, mini); }
  return 0;
}

/********* CONTAINER HOUSEKEEPING *********/

// Iterate over the raws in a bi and update reads/nraw
void bi_census(Bi *bi) {
  unsigned int r, reads=0, nraw=0;
  for(r=0;r<bi->nraw;r++) {
    reads += bi->raw[r]->reads;
    nraw++;
  }
  if(reads != bi->reads) {
    bi->update_e = true;
  }
  bi->reads = reads;
  bi->nraw = nraw;
}

// Takes a Bi object, and calculates and assigns its center Raw.
// Currently this is done trivially by choosing the most abundant raw.
// This function also currently assigns the cluster sequence (equal to center->seq).
void bi_assign_center(Bi *bi) {
  unsigned int r, max_reads;
  
  // Assign the raw with the most reads as the center
  bi->center = NULL;
  for(r=0,max_reads=0;r<bi->nraw;r++) {
    bi->raw[r]->lock = false; // Unlock everything as center is changing
    if(bi->raw[r]->reads > max_reads) { // Most abundant
      bi->center = bi->raw[r];
      max_reads = bi->center->reads;
    }
  }
  // Assign center sequence to bi->seq and flag check_locks
  if(bi->center) { strcpy(bi->seq, bi->center->seq); }
  bi->check_locks = true;
}

