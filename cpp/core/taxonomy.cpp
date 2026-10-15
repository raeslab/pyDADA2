// Adapted from DADA2 src/taxonomy.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * C_assign_taxonomy2 split into TaxonomyModel (reference k-mer model, fit once) and
//    assign_taxonomy (classification), using dada2::parallel_for.
//  * Bootstrap uniforms: Rcpp::runif (R's RNG) -> seeded std::mt19937_64, drawn up front in the
//    same count and order as upstream (nseq * NBOOT * (max_arraylen/8)). Converted to [0,1) with
//    (x >> 11) * 2^-53 so the stream is identical on every platform.
//  * Tie-breaking in get_best_genus: upstream seeded a fresh std::mt19937 from
//    std::random_device on every call (non-reproducible); now a per-sequence std::mt19937_64
//    seeded from (seed, sequence index), so results depend only on (inputs, seed).
//  * Fixed-size stack arrays (karray[9999], bootarray[9999/8]) -> std::vector sized per
//    sequence (upstream overflowed for sequences longer than 9999 nt).
//  * Outputs are 0-based with -1 for NA.

#include <algorithm>
#include <random>

#include "api.h"
#include "dada.h"
#include "parallel.h"

namespace {

// Gets kmer index
// Returns -1 if non-ACGT base encountered
int tax_kmer(const char *seq, unsigned int k) {
  unsigned int j, nti;
  int kmer=0;

  for(j=0; j<k; j++) {
    if(seq[j] == 'A') {
      nti = 0;
    } else if (seq[j] == 'C') {
      nti = 1;
    } else if (seq[j] == 'G') {
      nti = 2;
    } else if (seq[j] == 'T') {
      nti = 3;
    } else {
      kmer = -1;
      break;
    }
    kmer = 4*kmer + nti;
  }
  return(kmer);
}

// Sets to 1 (TRUE) the value of kvec corresponding to each valid kmer index in the provided sequence
void tax_kvec(const char *seq, unsigned int k, unsigned char *kvec) {
  unsigned int i;
  unsigned int len = strlen(seq);
  size_t klen = len - k + 1; // The number of kmers in this sequence
  int kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  for(i=0;i<n_kmers;i++) { kvec[i] = 0; }
  if(len < k) { return; } // pydada2: upstream's klen underflowed for references shorter than k

  for(i=0; i<klen; i++) {
    kmer = tax_kmer(&seq[i], k);

    // Ensure a valid kmer index
    if(kmer>=0 && (size_t)kmer<n_kmers) {
      kvec[kmer] = 1;
    }
  }
}

// Writes all valid (>=0) kmer indices in the provided sequence to karray. Returns number written.
unsigned int tax_karray(const char *seq, unsigned int k, int *karray) {
  unsigned int i, j;
  int kmer;
  unsigned int len = strlen(seq);
  if(len < k) { return 0; }
  size_t klen = len - k + 1; // The number of kmers in this sequence

  for(i=0,j=0;i<klen;i++) {
    kmer = tax_kmer(&seq[i], k);
    // Ensure a valid kmer index
    if(kmer>=0) {
      karray[j] = kmer;
      j++;
    }
  }
  std::sort(karray, karray+j);
  return(j);
}

// Uniform double in [0, 1) from a 64-bit generator, identical on all platforms.
inline double unif01(std::mt19937_64 &gen) {
  return (double)(gen() >> 11) * (1.0 / 9007199254740992.0);
}

int get_best_genus(const int *karray, float *out_logp, unsigned int arraylen, unsigned int n_kmers, unsigned int ngenus, const float *lgk_probability, std::mt19937_64 &gen) {
  unsigned int pos;
  const float *lgk_v;
  int kmer, max_g = -1;
  unsigned int g;
  float logp, max_logp = -FLT_MAX; // Init value to be replaced on first iteration
  double rv; // Dummy random variable
  unsigned int nmax=0; // Number of times the current max logp has been seen

  for(g=0;g<ngenus;g++) {
    lgk_v = &lgk_probability[(size_t)g*n_kmers];
    logp = 0.0;

    // Take the product of the probabilitys -> sum of logs
    // This is the rate limiting step of the entire assignTaxonomy (on query sets of non-trival size)
    for(pos=0;pos<arraylen;pos++) {
      kmer = karray[pos];
      logp += lgk_v[kmer];
      if(logp < max_logp) { break; }
    }

    if(max_logp > 0 || logp>max_logp) { // Store if new max
      max_logp = logp;
      max_g = g;
      nmax=1;
    } else if (max_logp == logp) { // With uniform prob, store if equal to current max
      nmax++;
      rv = unif01(gen);
      if(rv < 1.0/nmax) {
        max_g = g;
      }
    }
  }
  *out_logp = max_logp;
  return max_g;
}

// SplitMix64 finaliser: decorrelates (seed, index) pairs for per-sequence tie-break streams.
inline uint64_t mix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

constexpr unsigned int K = 8;

}  // namespace

namespace dada2 {

TaxonomyModel::TaxonomyModel(const std::vector<std::string> &refs, const std::vector<int> &ref_to_genus,
                             const std::vector<int> &genusmat, int ngenus, int nlevel)
    : ngenus_(ngenus), nlevel_(nlevel), genusmat_(genusmat) {
  size_t i, g;
  unsigned int k = K;
  size_t n_kmers = (1 << (2*k));
  size_t nref = refs.size();
  if(nref == 0) dada2_stop_invalid("No reference sequences provided.");
  if(nref != ref_to_genus.size()) dada2_stop_invalid("Length mismatch between number of references and map to genus.");
  if(ngenus < 1 || nlevel < 1) dada2_stop_invalid("Invalid genus matrix dimensions.");
  if(genusmat.size() != (size_t)ngenus * nlevel) dada2_stop_invalid("Genus matrix has the wrong size.");

  for(i=0;i<ref_to_genus.size();i++) {
    if(ref_to_genus[i]<0 || ref_to_genus[i] >= ngenus) {
      dada2_stop_invalid("Invalid map from references to genus.");
    }
  }

  // Count seqs in each genus (M_g).
  std::vector<float> genus_num_plus1(ngenus, 0.0f);
  for(i=0;i<nref;i++) {
    genus_num_plus1[ref_to_genus[i]]++;
  }
  for(g=0;g<(size_t)ngenus;g++) {
    genus_num_plus1[g]++;
  }

  std::vector<float> kmer_prior(n_kmers, 0.0f);
  lgk_probability_.assign((size_t)ngenus * n_kmers, 0.0f);
  float *lgk_v;
  std::vector<unsigned char> ref_kv(n_kmers);

  for(i=0;i<nref;i++) {
    // Calculate kmer-vector of this reference sequences
    tax_kvec(refs[i].c_str(), k, ref_kv.data());
    // Assign the kmer-counts to the appropriate "genus" and kmer-prior
    g = ref_to_genus[i];
    lgk_v = &lgk_probability_[g*n_kmers];
    for(size_t kmer=0;kmer<n_kmers;kmer++) {
      if(ref_kv[kmer]) {
        lgk_v[kmer]++;
        kmer_prior[kmer]++;
      }
    }
  }

  // Correct word priors
  for(size_t kmer=0;kmer<n_kmers;kmer++) {
    kmer_prior[kmer] = (kmer_prior[kmer] + 0.5)/(1.0 + nref);
  }

  ///! Create log genus-kmer probability
  for(g=0;g<(size_t)ngenus;g++) {
    lgk_v = &lgk_probability_[g*n_kmers];
    for(size_t kmer=0;kmer<n_kmers;kmer++) {
      lgk_v[kmer] = logf((lgk_v[kmer] + kmer_prior[kmer])/genus_num_plus1[g]);
    }
  }
}

TaxonomyResult assign_taxonomy(const TaxonomyModel &model, const std::vector<std::string> &seqs,
                               const std::vector<std::string> &rcs, bool try_rc, uint64_t seed,
                               int n_threads) {
  size_t i;
  unsigned int k = K;
  size_t n_kmers = (1 << (2*k));
  size_t nseq = seqs.size();
  if(nseq == 0) dada2_stop_invalid("No seqs provided to classify.");
  if(try_rc && rcs.size() != nseq) dada2_stop_invalid("Length mismatch between sequences and reverse complements.");
  const size_t ngenus = model.ngenus();
  const size_t nlevel = model.nlevel();
  const float *lgk_probability = model.log_probability().data();
  const int *C_genusmat = model.genusmat().data();

  // Get size of the kmer arrays for the sequences to be classified.
  unsigned int max_arraylen = 0;
  unsigned int seqlen;
  for(i=0;i<nseq;i++) {
    seqlen = seqs[i].size();
    if(seqlen < k) { continue; } // pydada2: upstream's unsigned wrap here made the allocation explode
    if((seqlen-k+1) > max_arraylen) { max_arraylen = seqlen-k+1; }
  }

  // Generate random numbers for bootstrapping.
  std::mt19937_64 unif_gen(seed);
  std::vector<double> C_unifs((size_t)nseq*NBOOT*(max_arraylen/8));
  for(i=0;i<C_unifs.size();i++) { C_unifs[i] = unif01(unif_gen); }

  TaxonomyResult res;
  res.tax.assign(nseq, -1);
  res.boot.assign(nseq * nlevel, 0);
  res.boot_tax.assign(nseq * NBOOT, -1);
  int *C_rval = res.tax.data();
  int *C_rboot = res.boot.data();
  int *C_rboot_tax = res.boot_tax.data();

  auto body = [&](std::size_t begin, std::size_t end) {
    size_t i, seqlen;
    unsigned int boot, booti, arraylen, arraylen_rc;
    int max_g, max_g_rc, boot_g;
    const double *unifs;
    float logp, logp_rc;
    std::vector<int> karray, karray_rc, bootarray;

    for(std::size_t j=begin;j<end;j++) {
      seqlen = seqs[j].size();
      if(seqlen < 50) { // No assignment made for very short seqeunces
        // NA assignments and 0 bootstrap confidences for this sequence (preset above)
        continue;
      }
      std::mt19937_64 tie_gen(mix64(seed ^ mix64((uint64_t)j)));
      karray.resize(seqlen);
      arraylen = tax_karray(seqs[j].c_str(), k, karray.data());

      // Find best hit
      max_g = get_best_genus(karray.data(), &logp, arraylen, n_kmers, ngenus, lgk_probability, tie_gen);
      if(try_rc) { // see if rev-comp is a better match to refs
        karray_rc.resize(rcs[j].size() > seqlen ? rcs[j].size() : seqlen);
        arraylen_rc = tax_karray(rcs[j].c_str(), k, karray_rc.data());
        if(arraylen != arraylen_rc) { dada2_stop("Discrepancy between forward and RC arraylen."); }
        max_g_rc = get_best_genus(karray_rc.data(), &logp_rc, arraylen_rc, n_kmers, ngenus, lgk_probability, tie_gen);
        if(logp_rc > logp) { // rev-comp is better, replace with it
          max_g = max_g_rc;
          memcpy(karray.data(), karray_rc.data(), arraylen * sizeof(int));
        }
      }

      C_rval[j] = max_g;

      unifs = &C_unifs[j*max_arraylen];
      booti = 0;
      bootarray.resize(arraylen/8 + 1);
      for(boot=0;boot<(unsigned int)NBOOT;boot++) {
        for(i=0;i<(arraylen/8);i++,booti++) {
          bootarray[i] = karray[(int) (arraylen*unifs[booti])];
        }
        boot_g = get_best_genus(bootarray.data(), &logp, (arraylen/8), n_kmers, ngenus, lgk_probability, tie_gen);
        C_rboot_tax[j*NBOOT+boot] = boot_g;
        for(i=0;i<nlevel;i++) {
          if(C_genusmat[boot_g*nlevel+i] == C_genusmat[max_g*nlevel+i]) {
            C_rboot[j*nlevel+i]++;
          } else {
            break;
          }
        }
      } // for(boot=0;boot<NBOOT;boot++)
    } // for(std::size_t j=begin;j<end;j++)
  };

  const size_t INTERRUPT_BLOCK_SIZE = 128;
  for(i=0;i<nseq;i+=INTERRUPT_BLOCK_SIZE) {
    size_t j = i+INTERRUPT_BLOCK_SIZE;
    if(j > nseq) { j = nseq; }
    dada2::parallel_for(i, j, body, 1, n_threads); // GRAIN_SIZE=1
    dada2_check_interrupt();
  }
  return res;
}

}  // namespace dada2
