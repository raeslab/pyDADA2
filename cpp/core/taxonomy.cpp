// Adapted from DADA2 src/taxonomy.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * C_assign_taxonomy2 split into TaxonomyModel (reference k-mer model, fit once) and
//    assign_taxonomy (classification), using dada2::parallel_for.
//  * Model build (same table values as upstream): upstream zeroed and scanned a 4^8-entry
//    presence vector per reference and incremented float counts per genus, serially. Here each
//    reference's distinct 8-mers are found with a rolling 2-bit code, counted as integers in
//    parallel over blocks of genera and converted to float (counts saturate at 2^24 exactly like
//    repeated float increments), then the same float expressions give the prior and the
//    log-probabilities.
//  * Table layout k-mer-major ([kmer][genus], genera padded to 16) instead of genus-major, and
//    get_best_genus sums the rows of the query's k-mers for all genera at once instead of one
//    genus at a time with an early exit. Each genus's sum is formed in the same order, and every
//    log-probability is < 0, so a partial sum that triggered upstream's early exit is >= the full
//    sum: the arg-max, the ties and the tie-breaking draws are identical.
//  * Bootstrap uniforms: Rcpp::runif (R's RNG) -> seeded std::mt19937_64, drawn up front in the
//    same order as upstream. Converted to [0,1) with (x >> 11) * 2^-53 so the stream is identical
//    on every platform. Upstream drew nseq * NBOOT * (max_arraylen/8) values but only reads
//    indices < (nseq-1) * max_arraylen + NBOOT * (max_arraylen/8); only that prefix is drawn.
//  * Tie-breaking in get_best_genus: upstream seeded a fresh std::mt19937 from
//    std::random_device on every call (non-reproducible); now a per-sequence std::mt19937_64
//    seeded from (seed, sequence index), so results depend only on (inputs, seed).
//  * Fixed-size stack arrays (karray[9999], bootarray[9999/8]) -> std::vector sized per
//    sequence (upstream overflowed for sequences longer than 9999 nt).
//  * Outputs are 0-based with -1 for NA.

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <mutex>
#include <new>
#include <random>

#include "api.h"
#include "dada.h"
#include "parallel.h"

namespace {

constexpr unsigned int K = 8;
constexpr size_t N_KMERS = size_t(1) << (2 * K);  // 4^8
constexpr size_t GENUS_BLOCK = 16;                // genera per block: 64 bytes of floats
constexpr uint32_t FLOAT_COUNT_MAX = 16777216;    // 2^24: float++ saturates here
constexpr std::align_val_t TABLE_ALIGN{64};

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

// Calls fn(kmer) for every valid 8-mer (all of A/C/G/T) of seq, left to right. Equivalent to
// tax_kmer at every position, with a rolling code.
template <typename Fn>
void for_each_kmer(std::string_view seq, Fn &&fn) {
  uint32_t code = 0;
  unsigned int valid = 0;  // length of the current run of A/C/G/T
  for(char c : seq) {
    uint32_t nti;
    switch(c) {
      case 'A': nti = 0; break;
      case 'C': nti = 1; break;
      case 'G': nti = 2; break;
      case 'T': nti = 3; break;
      default: valid = 0; continue;
    }
    code = ((code << 2) | nti) & (uint32_t)(N_KMERS - 1);
    if(++valid >= K) fn(code);
  }
}

// Uniform double in [0, 1) from a 64-bit generator, identical on all platforms.
inline double unif01(std::mt19937_64 &gen) {
  return (double)(gen() >> 11) * (1.0 / 9007199254740992.0);
}

// Upstream get_best_genus on the k-mer-major table: acc[g] = sum over pos of logp(karray[pos], g)
// in pos order (4 rows per pass over acc, same per-genus order), then upstream's sequential
// arg-max over g with its tie-breaking draws. acc has ngp entries.
int get_best_genus(const int *karray, float *out_logp, unsigned int arraylen, size_t ngp,
                   unsigned int ngenus, const float *table, float *acc, std::mt19937_64 &gen) {
  std::fill(acc, acc + ngp, 0.0f);
  unsigned int pos = 0;
  for(; pos + 4 <= arraylen; pos += 4) {
    const float *r0 = &table[(size_t)karray[pos] * ngp];
    const float *r1 = &table[(size_t)karray[pos+1] * ngp];
    const float *r2 = &table[(size_t)karray[pos+2] * ngp];
    const float *r3 = &table[(size_t)karray[pos+3] * ngp];
    for(size_t g=0; g<ngp; g++) {
      float a = acc[g];
      a += r0[g];
      a += r1[g];
      a += r2[g];
      a += r3[g];
      acc[g] = a;
    }
  }
  for(; pos < arraylen; pos++) {
    const float *r0 = &table[(size_t)karray[pos] * ngp];
    for(size_t g=0; g<ngp; g++) acc[g] += r0[g];
  }

  int max_g = -1;
  float logp, max_logp = -FLT_MAX; // Init value to be replaced on first iteration
  double rv; // Dummy random variable
  unsigned int nmax=0; // Number of times the current max logp has been seen
  for(unsigned int g=0;g<ngenus;g++) {
    logp = acc[g];
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

float *alloc_table(size_t n) {
  return static_cast<float *>(::operator new[](n * sizeof(float), TABLE_ALIGN));
}

}  // namespace

namespace dada2 {

void TaxonomyModel::AlignedDelete::operator()(float *p) const {
  ::operator delete[](p, TABLE_ALIGN);
}

size_t TaxonomyModel::table_size() const { return N_KMERS * (size_t)ngenus_pad_; }

TaxonomyModel::TaxonomyModel(const std::vector<std::string_view> &refs,
                             const std::vector<int> &ref_to_genus, const std::vector<int> &genusmat,
                             int ngenus, int nlevel, int n_threads)
    : ngenus_(ngenus), nlevel_(nlevel), ngenus_pad_(0), genusmat_(genusmat) {
  size_t i, g;
  const size_t nref = refs.size();
  if(nref == 0) dada2_stop_invalid("No reference sequences provided.");
  if(nref != ref_to_genus.size()) dada2_stop_invalid("Length mismatch between number of references and map to genus.");
  if(ngenus < 1 || nlevel < 1) dada2_stop_invalid("Invalid genus matrix dimensions.");
  if(genusmat.size() != (size_t)ngenus * nlevel) dada2_stop_invalid("Genus matrix has the wrong size.");

  for(i=0;i<ref_to_genus.size();i++) {
    if(ref_to_genus[i]<0 || ref_to_genus[i] >= ngenus) {
      dada2_stop_invalid("Invalid map from references to genus.");
    }
  }

  // Count seqs in each genus (M_g), as upstream (float increments).
  std::vector<float> genus_num_plus1(ngenus, 0.0f);
  for(i=0;i<nref;i++) {
    genus_num_plus1[ref_to_genus[i]]++;
  }
  for(g=0;g<(size_t)ngenus;g++) {
    genus_num_plus1[g]++;
  }

  // References grouped by genus (counting sort, stable).
  std::vector<size_t> goff((size_t)ngenus + 1, 0), gref(nref);
  for(i=0;i<nref;i++) goff[(size_t)ref_to_genus[i] + 1]++;
  for(g=0;g<(size_t)ngenus;g++) goff[g+1] += goff[g];
  {
    std::vector<size_t> pos(goff.begin(), goff.end() - 1);
    for(i=0;i<nref;i++) gref[pos[ref_to_genus[i]]++] = i;
  }

  const size_t ngp = ((size_t)ngenus + GENUS_BLOCK - 1) / GENUS_BLOCK * GENUS_BLOCK;
  ngenus_pad_ = (int)ngp;
  table_.reset(alloc_table(N_KMERS * ngp));
  float *T = table_.get();
  const size_t nblocks = ngp / GENUS_BLOCK;

  // Pass 1 (parallel over blocks of genera): number of references of each genus containing each
  // k-mer, written to the table as float, and the number of references containing each k-mer.
  std::vector<uint64_t> prior_count(N_KMERS, 0);
  std::mutex prior_mutex;
  auto count_body = [&](size_t b0, size_t b1) {
    std::vector<uint32_t> cnt(GENUS_BLOCK * N_KMERS);
    std::vector<uint32_t> prior(N_KMERS, 0);
    std::vector<uint32_t> stamp(N_KMERS, UINT32_MAX);  // last reference that contained a k-mer
    for(size_t blk=b0; blk<b1; blk++) {
      for(size_t gl=0; gl<GENUS_BLOCK; gl++) {
        const size_t gg = blk * GENUS_BLOCK + gl;
        if(gg >= (size_t)ngenus) break;
        uint32_t *c = &cnt[gl * N_KMERS];
        for(size_t r=goff[gg]; r<goff[gg+1]; r++) {
          const uint32_t ri = (uint32_t)gref[r];
          for_each_kmer(refs[ri], [&](uint32_t kmer) {
            if(stamp[kmer] != ri) {  // presence, not multiplicity (upstream tax_kvec)
              stamp[kmer] = ri;
              c[kmer]++;
              prior[kmer]++;
            }
          });
        }
      }
      for(size_t kmer=0; kmer<N_KMERS; kmer++) {  // write the block's columns, reset counts
        float *dst = &T[kmer * ngp + blk * GENUS_BLOCK];
        for(size_t gl=0; gl<GENUS_BLOCK; gl++) {
          uint32_t &cv = cnt[gl * N_KMERS + kmer];
          dst[gl] = (float)std::min(cv, FLOAT_COUNT_MAX);
          cv = 0;
        }
      }
    }
    std::lock_guard<std::mutex> lock(prior_mutex);
    for(size_t kmer=0; kmer<N_KMERS; kmer++) prior_count[kmer] += prior[kmer];
  };
  dada2::parallel_for(0, nblocks, count_body, 1, n_threads);

  // Correct word priors (upstream: float counts, then this expression)
  std::vector<float> kmer_prior(N_KMERS);
  for(size_t kmer=0;kmer<N_KMERS;kmer++) {
    kmer_prior[kmer] = (float)std::min<uint64_t>(prior_count[kmer], FLOAT_COUNT_MAX);
    kmer_prior[kmer] = (kmer_prior[kmer] + 0.5)/(1.0 + nref);
  }

  // Pass 2 (parallel over k-mers): log genus-kmer probability, as upstream.
  auto log_body = [&](size_t k0, size_t k1) {
    for(size_t kmer=k0; kmer<k1; kmer++) {
      float *row = &T[kmer * ngp];
      const float kp = kmer_prior[kmer];
#if defined(_MSC_VER) && !defined(__clang__)
#pragma loop(no_vector)  // keep the scalar CRT logf (MSVC may otherwise use a vector logf)
#endif
      for(size_t gg=0; gg<(size_t)ngenus; gg++) {
        row[gg] = logf((row[gg] + kp)/genus_num_plus1[gg]);
      }
      for(size_t gg=(size_t)ngenus; gg<ngp; gg++) row[gg] = 0.0f;
    }
  };
  dada2::parallel_for(0, N_KMERS, log_body, 256, n_threads);
}

TaxonomyModel::TaxonomyModel(const std::vector<int> &genusmat, int ngenus, int nlevel,
                             const float *table, size_t table_size)
    : ngenus_(ngenus), nlevel_(nlevel), ngenus_pad_(0), genusmat_(genusmat) {
  if(ngenus < 1 || nlevel < 1) dada2_stop_invalid("Invalid genus matrix dimensions.");
  if(genusmat.size() != (size_t)ngenus * nlevel) dada2_stop_invalid("Genus matrix has the wrong size.");
  const size_t ngp = ((size_t)ngenus + GENUS_BLOCK - 1) / GENUS_BLOCK * GENUS_BLOCK;
  if(table == nullptr || table_size != N_KMERS * ngp) dada2_stop_invalid("Taxonomy table has the wrong size.");
  ngenus_pad_ = (int)ngp;
  table_.reset(alloc_table(table_size));
  std::memcpy(table_.get(), table, table_size * sizeof(float));
}

TaxonomyResult assign_taxonomy(const TaxonomyModel &model, const std::vector<std::string> &seqs,
                               const std::vector<std::string> &rcs, bool try_rc, uint64_t seed,
                               int n_threads) {
  size_t i;
  unsigned int k = K;
  size_t nseq = seqs.size();
  if(nseq == 0) dada2_stop_invalid("No seqs provided to classify.");
  if(try_rc && rcs.size() != nseq) dada2_stop_invalid("Length mismatch between sequences and reverse complements.");
  const size_t ngenus = model.ngenus();
  const size_t nlevel = model.nlevel();
  const size_t ngp = model.ngenus_pad();
  const float *table = model.log_probability();
  const int *C_genusmat = model.genusmat().data();

  // Get size of the kmer arrays for the sequences to be classified.
  unsigned int max_arraylen = 0;
  unsigned int seqlen;
  for(i=0;i<nseq;i++) {
    seqlen = seqs[i].size();
    if(seqlen < k) { continue; } // pydada2: upstream's unsigned wrap here made the allocation explode
    if((seqlen-k+1) > max_arraylen) { max_arraylen = seqlen-k+1; }
  }

  // Generate random numbers for bootstrapping: the prefix of upstream's stream that is read
  // (sequence j reads indices j*max_arraylen + [0, NBOOT*(arraylen_j/8))).
  const size_t n_unifs_upstream = (size_t)nseq*NBOOT*(max_arraylen/8);
  const size_t n_unifs_read = (nseq-1)*(size_t)max_arraylen + (size_t)NBOOT*(max_arraylen/8);
  std::mt19937_64 unif_gen(seed);
  std::vector<double> C_unifs(std::min(n_unifs_upstream, n_unifs_read));
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
    std::vector<float> acc(ngp);

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
      max_g = get_best_genus(karray.data(), &logp, arraylen, ngp, ngenus, table, acc.data(), tie_gen);
      if(try_rc) { // see if rev-comp is a better match to refs
        karray_rc.resize(rcs[j].size() > seqlen ? rcs[j].size() : seqlen);
        arraylen_rc = tax_karray(rcs[j].c_str(), k, karray_rc.data());
        if(arraylen != arraylen_rc) { dada2_stop("Discrepancy between forward and RC arraylen."); }
        max_g_rc = get_best_genus(karray_rc.data(), &logp_rc, arraylen_rc, ngp, ngenus, table, acc.data(), tie_gen);
        if(logp_rc > logp) { // rev-comp is better, replace with it
          max_g = max_g_rc;
          memcpy(karray.data(), karray_rc.data(), arraylen * sizeof(int));
        }
      }

      C_rval[j] = max_g;

      unifs = (arraylen/8) > 0 ? C_unifs.data() + j*max_arraylen : nullptr;
      booti = 0;
      bootarray.resize(arraylen/8 + 1);
      for(boot=0;boot<(unsigned int)NBOOT;boot++) {
        for(i=0;i<(arraylen/8);i++,booti++) {
          bootarray[i] = karray[(int) (arraylen*unifs[booti])];
        }
        boot_g = get_best_genus(bootarray.data(), &logp, (arraylen/8), ngp, ngenus, table, acc.data(), tie_gen);
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
