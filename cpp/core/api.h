// Public C++ API of the pydada2 core (used by the nanobind bindings and C++ tests).
// One function per upstream Rcpp export (see DADA2 src/RcppExports.cpp), with plain C++ types.
//
// Conventions:
//  * Matrices are passed as contiguous ROW-major buffers with explicit dimensions.
//  * All indices are 0-based; "NA" is represented by -1 for integers and NaN for doubles.
//  * Invalid input -> std::invalid_argument; internal failures -> dada2::Error.

#ifndef DADA2_API_H
#define DADA2_API_H

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace dada2 {

// ------------------------------------------------------------------ dada_uniques (Rmain.cpp) ---
struct DadaOptions {
  int match = 5;
  int mismatch = -4;
  int gap = -8;
  bool use_kmers = true;
  double kdist_cutoff = 0.42;
  int band_size = 16;
  double omega_a = 1e-40;
  double omega_p = 1e-4;
  double omega_c = 1e-40;
  bool detect_singletons = false;
  int max_clust = 0;
  double min_fold = 1;
  int min_hamming = 1;
  int min_abund = 1;
  bool use_quals = true;
  bool final_consensus = false;
  bool vectorized_alignment = true;
  int homo_gap = -8;
  int sse = 2;
  bool gapless = true;
  bool greedy = true;
};

struct DadaOutput {
  // $clustering (one entry per cluster)
  std::vector<std::string> sequence;
  std::vector<int> abundance, n0, n1, nunq;
  std::vector<double> pval;
  std::vector<int> birth_from;          // 0-based, -1 for the first cluster
  std::vector<double> birth_pval, birth_fold;
  std::vector<int> birth_ham;           // -1 for the first cluster
  std::vector<double> birth_qave;       // NaN when not available
  // $birth_subs (one entry per birth substitution)
  std::vector<int> bs_pos;              // 0-based position in the parent center
  std::vector<std::string> bs_ref, bs_sub;
  std::vector<double> bs_qual;          // NaN without qualities
  std::vector<int> bs_clust;            // 0-based cluster index
  // $subqual: 16 x trans_ncol, row-major
  std::vector<int> trans;
  int trans_ncol = 0;
  // $clusterquals: nclust x quals_ncol (row-major; NaN beyond each center's length)
  std::vector<double> cluster_quals;
  int quals_ncol = 0;
  // $map: unique index -> cluster index (0-based), -1 when not corrected
  std::vector<int> map;
  // $pval: per-unique abundance p-value
  std::vector<double> pval_raw;
};

// err: 16 x err_ncol row-major. quals: nraw x quals_ncol row-major (quals_ncol == max sequence
// length), or nullptr/quals_ncol == 0 when there are no qualities.
DadaOutput dada_uniques(const std::vector<std::string> &seqs, const std::vector<int> &abundances,
                        const std::vector<bool> &priors, const double *err, int err_ncol,
                        const double *quals, int quals_ncol, const DadaOptions &opt, int n_threads,
                        bool verbose);

// ------------------------------------------------------------------ evaluate.cpp ---
std::pair<std::string, std::string> nwalign(const std::string &s1, const std::string &s2, int match,
                                            int mismatch, int gap_p, int homo_gap_p, int band,
                                            bool endsfree);
std::vector<std::pair<std::string, std::string>> nwvec(const std::vector<std::string> &s1,
                                                       const std::vector<std::string> &s2,
                                                       int match, int mismatch, int gap_p,
                                                       int band, bool endsfree);
std::array<int, 3> eval_pair(const std::string &s1, const std::string &s2);
std::string pair_consensus(const std::string &s1, const std::string &s2, int prefer,
                           bool trim_overhang);
std::vector<bool> is_acgt(const std::vector<std::string> &seqs);
std::vector<double> kmer_dist(const std::vector<std::string> &s1, const std::vector<std::string> &s2,
                              int kmer_size);
std::vector<double> kord_dist(const std::vector<std::string> &s1, const std::vector<std::string> &s2,
                              int kmer_size, int sse);
std::vector<int> kmer_matches(const std::vector<std::string> &s1,
                              const std::vector<std::string> &s2, int kmer_size);
std::vector<int> kdist_matches(const std::vector<std::string> &s1,
                               const std::vector<std::string> &s2, int kmer_size);

// ------------------------------------------------------------------ filter.cpp ---
std::vector<int> match_ref(const std::vector<std::string> &seqs, const std::string &ref,
                           int word_size, bool non_overlapping);
// Quality matrix nrow x ncol row-major; values == NA_QUAL terminate a row (like R's NA).
constexpr int32_t NA_QUAL = INT32_MIN;
std::vector<double> matrix_ee(const int32_t *quals, int nrow, int ncol);

// ------------------------------------------------------------------ chimera.cpp ---
bool is_bimera(const std::string &sq, const std::vector<std::string> &pars, bool allow_one_off,
               int min_one_off_par_dist, int match, int mismatch, int gap_p, int max_shift);
// mat: nsample x nseq row-major abundance table. Returns (nflag, nsam) per sequence.
std::pair<std::vector<int>, std::vector<int>> table_bimera(
    const int32_t *mat, int nrow, int ncol, const std::vector<std::string> &seqs, double min_fold,
    int min_abund, bool allow_one_off, int min_one_off_par_dist, int match, int mismatch,
    int gap_p, int max_shift, int n_threads);

// ------------------------------------------------------------------ taxonomy.cpp ---
constexpr int NBOOT = 100;

// Fit-once reference k-mer model (the first half of upstream C_assign_taxonomy2).
class TaxonomyModel {
 public:
  // ref_to_genus: 0-based genus index for each reference.
  // genusmat: ngenus x nlevel row-major integer taxonomy codes.
  TaxonomyModel(const std::vector<std::string> &refs, const std::vector<int> &ref_to_genus,
                const std::vector<int> &genusmat, int ngenus, int nlevel);
  int ngenus() const { return ngenus_; }
  int nlevel() const { return nlevel_; }
  const std::vector<float> &log_probability() const { return lgk_probability_; }
  const std::vector<int> &genusmat() const { return genusmat_; }

 private:
  int ngenus_, nlevel_;
  std::vector<float> lgk_probability_;
  std::vector<int> genusmat_;
};

struct TaxonomyResult {
  std::vector<int> tax;        // best genus per sequence (0-based), -1 when not assigned
  std::vector<int> boot;       // nseq x nlevel row-major bootstrap counts
  std::vector<int> boot_tax;   // nseq x NBOOT row-major bootstrap genus (0-based), -1 = NA
};

// The bootstrap uniforms and the tie-breaking stream are drawn from seeded generators, so the
// result depends only on (inputs, seed) and not on n_threads.
TaxonomyResult assign_taxonomy(const TaxonomyModel &model, const std::vector<std::string> &seqs,
                               const std::vector<std::string> &rcs, bool try_rc, uint64_t seed,
                               int n_threads);

}  // namespace dada2

#endif  // DADA2_API_H
