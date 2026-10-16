// C++ kernels ported from DADA2's R-level code (K1-K9, docs/porting-notes.md). Each kernel reproduces the R behaviour
// (including the Bioconductor helper semantics) of the function(s) it is ported from.

#ifndef DADA2_PORTED_H
#define DADA2_PORTED_H

#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dada2 {

// ------------------------------------------------------------------- K1 fastq_filter ---
// Per-mate filtering/trimming parameters (R: fastqFilter / fastqPairedFilter arguments).
struct MateFilterParams {
  int trunc_q = 2;           // truncate before the first base with numeric quality <= trunc_q
  int trunc_len = 0;         // 0 = no truncation
  double max_len = std::numeric_limits<double>::infinity();  // on the raw reads
  int min_len = 20;
  int trim_left = 0;
  int trim_right = 0;
  int max_n = 0;             // max number of non-ACGT letters
  int min_q = 0;             // applied only when min_q > trunc_q (as in R)
  double max_ee = std::numeric_limits<double>::infinity();
  bool rm_phix = true;
  double rm_lowcomplex = 0;  // 0 = off
  int qual_offset = 33;      // 33 (phred+33) or 64 (Solexa/phred+64)
};

struct MateFilterResult {
  std::vector<uint8_t> pass;  // per read
  std::vector<int> start;     // 0-based start of the kept subsequence
  std::vector<int> end;       // 0-based exclusive end
};

// Applies the per-mate pipeline (maxLen, trimLeft, trimRight, truncQ, truncLen, minLen, maxN,
// minQ, maxEE, phiX, low complexity) to already-oriented, upper-cased reads.
MateFilterResult filter_mate(const std::vector<std::string> &seqs,
                             const std::vector<std::string> &quals, const MateFilterParams &p,
                             const std::string &phix, int n_threads);

// ShortRead's .qualityTypeAuto on the first 10000 quality strings of a chunk:
// 64 (SFastqQuality) if the smallest char is > ':' and the largest > 'J', else 33.
int detect_quality_offset(const std::vector<std::string> &quals);

// isPhiX: word matches against the phiX genome and its reverse complement.
std::vector<bool> is_phix(const std::vector<std::string> &seqs, const std::string &phix,
                          int word_size, int min_matches, bool non_overlapping, int n_threads = 1);

// ------------------------------------------------------------------- K2 derep ---
// Streaming dereplication reproducing derepFastq/qtables2, including the chunk-dependent order
// of tied abundances (new uniques of each chunk appended in Biostrings lexical order, then a
// stable sort by decreasing abundance).
class Dereplicator {
 public:
  // Adds one chunk (as yielded by FastqStreamer with n records). Qualities are encoded strings.
  void add_chunk(const std::vector<std::string> &seqs, const std::vector<std::string> &quals,
                 int qual_offset);
  struct Result {
    std::vector<std::string> uniques;
    std::vector<int64_t> abundance;
    std::vector<double> quals;  // n_uniques x maxlen, row-major, NaN padded
    int maxlen = 0;
    std::vector<int64_t> map;   // read -> unique (0-based), -1 for zero-length reads
    int64_t n_zero_length = 0;
  };
  Result finish() const;

 private:
  std::vector<std::string> uniques_;
  std::vector<int64_t> counts_;
  std::vector<std::vector<double>> qsums_;
  std::vector<int64_t> map_;
  std::unordered_map<std::string, int64_t> index_;
  int maxlen_ = 0;
  int64_t n_zero_ = 0;
};

// Biostrings/ShortRead lexical order of DNA strings (IUPAC bit codes: A<C<M<G<R<S<V<T<W<Y<H<K<D<B<N).
bool dna_less(const std::string &a, const std::string &b);

// ------------------------------------------------------------------- K6 complexity ---
// seqComplexity: effective number of k-mers, exp(Shannon entropy) of the ACGT k-mer counts.
// window <= 0: whole sequence; otherwise the minimum over windows starting every `by` bases
// (initialised to 4^k, as in R). Sequences without any valid k-mer give NaN (R: NA).
std::vector<double> seq_complexity(const std::vector<std::string> &seqs, int kmer_size, int window,
                                   int by, int n_threads);

// ------------------------------------------------------------------- K7 primers ---
struct PrimerHits {
  std::vector<int> nhits;        // number of matches per sequence
  std::vector<int> first_start;  // 0-based start of the first match (negative: out of limits)
  std::vector<int> first_end;    // 0-based exclusive end of the first match
  std::vector<int> last_start;   // last match (largest start)
  std::vector<int> last_end;
};

// Biostrings matchPattern/vmatchPattern semantics: up to max_mismatch mismatches (positions
// outside the subject count as mismatches) or, with_indels, "best local matches" with edit
// distance <= max_mismatch. fixed=false: IUPAC codes in pattern and subject match when their
// base sets intersect; fixed=true: letters must be identical.
PrimerHits match_primer(const std::vector<std::string> &seqs, const std::string &primer,
                        int max_mismatch, bool with_indels, bool fixed, int n_threads);

// All matches (0-based start, exclusive end) of a pattern in one subject, sorted by start.
std::vector<std::pair<int, int>> match_pattern(const std::string &subject,
                                               const std::string &pattern, int max_mismatch,
                                               bool with_indels, bool fixed);

// ------------------------------------------------------------------- K3 merge ---
struct MergeResult {
  std::vector<int> nmatch, nmismatch, nindel;
  std::vector<std::string> sequence;  // pair consensus (gaps removed)
};

// Align fwd[i] with rev_rc[i] (nwalign, band=-1, ends-free), evaluate and build the consensus
// (prefer[i] = 1 or 2 wins mismatches).
MergeResult merge_pairs_core(const std::vector<std::string> &fwd, const std::vector<std::string> &rev_rc,
                             const std::vector<int> &prefer, int match, int mismatch, int gap,
                             bool trim_overhang, int n_threads);

// ------------------------------------------------------------------- K5 collapse / shift ---
// collapseNoMismatch's greedy loop over sequences in processing order; returns the index of the
// kept sequence each one is merged into (itself when kept).
std::vector<int> collapse_no_mismatch(const std::vector<std::string> &seqs, int min_overlap,
                                      bool vec, int band, int n_threads);

// isShiftDenovo: whether each sequence is a shifted copy of a more abundant one.
std::vector<bool> is_shift_denovo(const std::vector<std::string> &seqs,
                                  const std::vector<int64_t> &abunds, int min_overlap,
                                  bool flag_subseqs, int n_threads);

// ------------------------------------------------------------------- K4 bimera de novo ---
// isBimeraDenovo: for each unique, parents are the uniques with abundance > min_fold * own and
// > min_parent_abund (in input order); flagged if is_bimera() with >= 2 parents.
std::vector<bool> bimera_denovo(const std::vector<std::string> &seqs,
                                const std::vector<double> &abunds, double min_fold,
                                double min_parent_abund, bool allow_one_off,
                                int min_one_off_par_dist, int match, int mismatch, int gap_p,
                                int max_shift, int n_threads);

// ------------------------------------------------------------------- K8 species ---
// For each query, the (sorted, 0-based) indices of the references containing it exactly
// (also in reverse complement when try_rc).
// refs are views: the caller keeps the reference strings alive during the call.
std::vector<std::vector<int>> species_hits(const std::vector<std::string> &queries,
                                           const std::vector<std::string_view> &refs,
                                           bool try_rc, int n_threads);

}  // namespace dada2

#endif  // DADA2_PORTED_H
