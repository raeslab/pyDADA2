// Adapted from DADA2 src/error.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * Builders fill dada2::DadaOutput instead of returning Rcpp DataFrames/matrices.
//  * Indices are 0-based; R NA -> -1 (integers) / NaN (doubles).
//  * The cluster quality matrix is returned clusters x positions (R transposed it afterwards).
//  * Unused b_make_positional_substitution_df removed.
//  * Transition matrix: explicit range check on the quality index (upstream wrote out of bounds).
#include "output.h"


// This function constructs the output "clustering" data.frame for the dada(...) function.
// This contains core and diagnostic information on each partition (or cluster, or Bi).
void b_make_clustering_df(B *b, Sub **subs, Sub **birth_subs, bool has_quals, dada2::DadaOutput &out) {
  unsigned int i, j, r, s, cind, max_reads;
  Raw *max_raw, *raw;
  Sub *sub;
  double q_ave;
  size_t index;
  
  // Create output character-vector of representative sequences for each partition (Bi)
  std::vector<std::string> &Rseqs = out.sequence;
  Rseqs.clear();
  char oseq[SEQLEN];
  for(i=0;i<b->nclust;i++) {
    max_reads=0;
    max_raw = NULL;
    for(r=0;r<b->bi[i]->nraw;r++) {
      if(b->bi[i]->raw[r]->reads > max_reads) {
        max_raw = b->bi[i]->raw[r];
        max_reads = max_raw->reads;
      }
    }
    if(!max_raw) {
      Rseqs.push_back(std::string(""));
    } else {
      ntcpy(oseq, max_raw->seq);
      Rseqs.push_back(std::string(oseq));
    }
  }
  
  // Create output vectors for other columns of clustering data.frame
  // Each has one entry for each partition (Bi)
  auto &Rabunds = out.abundance;      // abundances
  Rabunds.assign(b->nclust, 0);
  auto &Rzeros = out.n0;       // n0
  Rzeros.assign(b->nclust, 0);
  auto &Rones = out.n1;        // n1
  Rones.assign(b->nclust, 0);
  auto &Rraws = out.nunq;        // nraw
  Rraws.assign(b->nclust, 0);
  auto &Rbirth_pvals = out.birth_pval; // pvalue at birth
  Rbirth_pvals.assign(b->nclust, 0);
  auto &Rbirth_folds = out.birth_fold; // fold over-abundance at birth
  Rbirth_folds.assign(b->nclust, 0);
  auto &Rbirth_hams = out.birth_ham;  // hamming distance at birth
  Rbirth_hams.assign(b->nclust, 0);
  std::vector<double> Rbirth_es(b->nclust);    // expected number at birth (not returned upstream either)
///  Rcpp::CharacterVector Rbirth_types;           // DEPRECATED
  auto &Rbirth_froms = out.birth_from; // cluster from which this new cluster was born
  Rbirth_froms.assign(b->nclust, 0);
  auto &Rbirth_qaves = out.birth_qave; // average quality of substitutions that drove birth
  Rbirth_qaves.assign(b->nclust, 0);
  auto &Rpvals = out.pval;       // post-hoc pvalue
  Rpvals.assign(b->nclust, 0);

  // Assign values to the output vectors
  for(i=0;i<b->nclust;i++) {
    // abundance, nraw, n0, n1
    Rzeros[i] = 0; Rones[i] = 0;
    for(r=0;r<b->bi[i]->nraw;r++) {
      raw = b->bi[i]->raw[r];
      if(raw->correct) { // Use only correct=true raws (can assume true of center, but 2x check)
        Rabunds[i] += raw->reads;
        Rraws[i]++;
        sub = subs[raw->index];
        if(sub) {
          if(sub->nsubs == 0) { Rzeros[i] += b->bi[i]->raw[r]->reads; }
          if(sub->nsubs == 1) { Rones[i] += b->bi[i]->raw[r]->reads; }
        }
      }
    }
    // Record information from the cluster's birth
///    Rbirth_types.push_back(std::string(b->bi[i]->birth_type));
    if(i==0) {  // 0-clust wasn't born normally
      Rbirth_pvals[i] = NAN;
      Rbirth_froms[i] = -1;
      Rbirth_folds[i] = NAN;
      Rbirth_hams[i] = -1;
      Rbirth_es[i] = NAN;
      Rbirth_qaves[i] = NAN;
    } else { 
      Rbirth_froms[i] = b->bi[i]->birth_from; // 0-based (upstream: +1 for R)
      Rbirth_pvals[i] = b->bi[i]->birth_pval;
      Rbirth_folds[i] = b->bi[i]->birth_fold;
      Rbirth_hams[i] = b->bi[i]->birth_comp.hamming;
      Rbirth_es[i] = b->bi[i]->birth_e;
      // Calculate average quality of birth substitutions
      if(has_quals) {
        q_ave = 0.0;
        sub = birth_subs[i];
        if(sub && sub->q1) {
          for(s=0;s<sub->nsubs;s++) {
            q_ave += (sub->q1[s]); // Numeric vector for output
          }
          q_ave = q_ave/((double)sub->nsubs);
        }
        Rbirth_qaves[i] = q_ave;
      } else {
        Rbirth_qaves[i] = NAN;
      }
    }
  }
  
  // Calculate post-hoc pval
  // This is not as exhaustive anymore
  std::unordered_map<unsigned int, unsigned int> center_of;
  for(i=0;i<b->nclust;i++) {
    center_of[b->bi[i]->center->index] = i;
  }
  std::vector<double> tot_e(b->nclust);
  for(i=0;i<b->nclust;i++) {
    for(cind=0;cind<b->bi[i]->comp.size();cind++) {
      index = b->bi[i]->comp[cind].index;
      if(center_of.count(index)) { // This is the center of a Bi
        j = center_of[index];
        if(i != j) { // Only add from other Bis
          tot_e[j] += b->bi[i]->comp[cind].lambda * b->bi[i]->reads;
        }
      }
    }
  }
  for(i=0;i<b->nclust;i++) {
    Rpvals[i] = calc_pA(b->bi[i]->center->reads, tot_e[i], true); // prior=true to get non-conditional p-val
  }
  
}

// Returns a 16xN matrix with the observed counts of each transition categorized by
// type (row) and quality (column). Assumes qualities start at 0.
void b_make_transition_by_quality_matrix(B *b, Sub **subs, bool has_quals, int ncol, dada2::DadaOutput &out) {
  unsigned int i, r, pos0, pos1, nti0, nti1, qual, t_ij;
  Sub *sub;
  Raw *raw, *center;
  
  if(!has_quals) { ncol = 1; }
  
  // Storage for counts for 0...(ncol-1) and each nti->ntj
  // 16 x ncol, row-major
  out.trans.assign(16 * (size_t)ncol, 0);
  out.trans_ncol = ncol;
  auto transMat = [&](unsigned int t, unsigned int q) -> int & { return out.trans[(size_t)t * ncol + q]; };

  for(i=0;i<b->nclust;i++) {
    center = b->bi[i]->center;
    for(r=0;r<b->bi[i]->nraw;r++) {
      raw = b->bi[i]->raw[r];
      if(!raw->correct) { continue; } // Don't count if not being corrected... Need to make sure used w/ care R side
      sub = subs[raw->index]; // The sub object includes the map between the center and the raw positions
      if(!sub) {
        if(VERBOSE) { dada2_log("Warning: No sub for R%i in C%i.\n", r, i); }
        continue;
      }
      
      for(pos0=0;pos0<center->length;pos0++) {
        pos1 = sub->map[pos0];
        if(pos1 == GAP_GLYPH) { // A gap in the aligned seq
          continue; // Gaps excluded from the model
        }
        nti0 = (int) (center->seq[pos0] - 1);
        nti1 = (int) (raw->seq[pos1] - 1);
        qual = raw->qual[pos1]; // qual=unsigned int
        // And record these counts
        t_ij = (4*nti0)+nti1;
        if(has_quals) {
          if((int)qual >= ncol) { dada2_stop("Rounded quality exceeded range of err lookup table."); }
          transMat(t_ij, qual) += raw->reads;
        } else { 
          transMat(t_ij, 0) += raw->reads; 
        }
      } // for(pos0=0;pos0<center->length;pos0++)
    } // for(r=0;b->bi[i]->nraw)
  } // for(i=0;i<b->nclust;i++)
  
}

// Calculate the average positional qualities for each cluster/partition/Bi
// Return Bi (rows) by position (columns) matrix, row-major (upstream: position x Bi, transposed in R).
void b_make_cluster_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int maxlen, dada2::DadaOutput &out) {
  unsigned int i, r, pos0, pos1, raw_reads, seqlen;
  std::vector<unsigned int> nreads(maxlen);
  Sub *sub;
  Raw *raw;
  out.cluster_quals.assign((size_t)maxlen * b->nclust, 0.0);
  out.quals_ncol = maxlen;
  auto Rquals = [&](unsigned int pos, unsigned int i) -> double & { return out.cluster_quals[(size_t)i * maxlen + pos]; };
  
  if(has_quals) {
    for(i=0;i<b->nclust;i++) {
      seqlen = b->bi[i]->center->length;
      for(pos0=0;pos0<seqlen;pos0++) { nreads[pos0] = 0; }
      for(r=0;r<b->bi[i]->nraw;r++) {
        raw = b->bi[i]->raw[r];
        if(!raw->correct) { continue; }
        raw_reads = raw->reads;
        sub = subs[raw->index];
        if(sub) {
          for(pos0=0;pos0<seqlen;pos0++) {
            pos1 = sub->map[pos0];
            if(pos1 == GAP_GLYPH) { // Gap
              continue;
            }
            nreads[pos0] += raw_reads;
            Rquals(pos0,i) += (raw->qual[pos1] * raw_reads);  // Output qual construction
          }
        }
      } // for(pos0=0;pos0<len1;pos0++)
      for(pos0=0;pos0<seqlen;pos0++) { Rquals(pos0,i) = Rquals(pos0,i)/nreads[pos0]; }
      for(pos0=seqlen;pos0<maxlen;pos0++) { Rquals(pos0,i) = NAN; }
    }
  }
}

// Make output data.frame of birth_subs
void b_make_birth_subs_df(B *b, Sub **birth_subs, bool has_quals, dada2::DadaOutput &out) {
  Sub *sub;
  unsigned int i, j, s;
  unsigned int tot_subs=0;
  char buf[2] = {'\0','\0'};
  // Count total number of subs
  for(i=0;i<b->nclust;i++) {
    sub = birth_subs[i];
    if(sub) { tot_subs += sub->nsubs; }
  }
  // Initialize the columns
  auto &bs_pos = out.bs_pos; bs_pos.assign(tot_subs, 0);
  auto &bs_nt0 = out.bs_ref; bs_nt0.assign(tot_subs, std::string());
  auto &bs_nt1 = out.bs_sub; bs_nt1.assign(tot_subs, std::string());
  auto &bs_qual = out.bs_qual; bs_qual.assign(tot_subs, 0.0);
  auto &bs_clust = out.bs_clust; bs_clust.assign(tot_subs, 0);
  // Record data for each birth substitution
  for(i=0,j=0;i<b->nclust;i++) {
    sub = birth_subs[i];
    if(sub) {
      for(s=0;s<sub->nsubs;s++) {
        bs_pos[j] = sub->pos[s]; // 0-based (upstream: +1 for R)
        buf[0] = sub->nt0[s];
        int2nt(buf, buf);
        bs_nt0[j].assign(std::string(buf));
        buf[0] = sub->nt1[s];
        int2nt(buf, buf);
        bs_nt1[j].assign(std::string(buf));
        if(has_quals) {
          bs_qual[j] = sub->q1[s];  // NumericVector, making birth sub quals from sub object
        } else {
          bs_qual[j] = NAN;
        }
        bs_clust[j] = i; // 0-based (upstream: +1 for R)
        j++;
      }
    }
  }
}
