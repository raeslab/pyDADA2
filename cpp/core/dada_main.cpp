// Adapted from DADA2 src/Rmain.cpp (https://github.com/benjjneb/dada2, commit b61fe2e, v1.36.0).
// Copyright (C) Benjamin Callahan et al. Licensed under the LGPL (see LICENSE).
//
// Changes for pydada2:
//  * dada_uniques takes/returns plain C++ types (api.h); 0-based indices, -1 for NA.
//  * multithread flag -> n_threads (dada2::parallel_for); b_compare/b_compare_parallel merged.
//  * Quality matrix is nraw x maxlen row-major (R passed the transpose, column-major).
//  * Additional input validation (non-ACGT, abundances, NaN/out-of-range qualities) that upstream
//    performed R-side.
//  * Exception-safe ownership of raws/B/subs (DadaState); Rcpp::checkUserInterrupt ->
//    dada2_check_interrupt (handler installed by the bindings).
//  * SSE is no longer forced to 0 on non-x86 (exact scalar emulation in kmers.cpp).
//  * Validation of the error matrix (finite, in [0, 1], as R's getErrors(enforce=TRUE)) and of
//    the presence of qualities when use_quals is set (upstream read a NULL quality array).
//  * The B object is owned by DadaState before the algorithm runs (it leaked when an exception,
//    e.g. an interrupt, escaped run_dada).
//  * Raw construction and the k-mer passes run in parallel over raws; the final substitutions
//    are computed in parallel over raws instead of clusters (cluster 0 holds most raws). Each
//    raw writes only its own slots, so the results do not depend on n_threads.
//  * Optional comparison cache across calls (AlignmentCache) and dada_uniques_batch, which runs
//    several samples concurrently on the shared thread pool.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>
#include <tuple>
#include <unordered_map>

#include "api.h"
#include "dada.h"
#include "output.h"
#include "parallel.h"

void run_dada(B *bb, const ErrMat &errMat,
              int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size,
              bool detect_singletons, int max_clust, double min_fold, int min_hamming, int min_abund,
              bool vectorized_alignment, int n_threads, bool verbose, int SSE, bool gapless, bool greedy);

namespace {

// Owns the C-style objects of one dada_uniques call so that they are released when an
// exception propagates (upstream relied on R's longjmp and leaked on Rcpp::stop).
struct DadaState {
  std::vector<Raw *> raws;
  B *bb = nullptr;
  std::vector<Sub *> subs;
  std::vector<Sub *> birth_subs;
  std::vector<uint8_t> k8;
  std::vector<uint16_t> k16;
  std::vector<uint16_t> kord;
  ~DadaState() {
    for(Sub *s : subs) sub_free(s);
    for(Sub *s : birth_subs) sub_free(s);
    if(bb) b_free(bb);
    for(Raw *r : raws) if(r) raw_free(r);
  }
};

}  // namespace

/********** COMPARISON CACHE *********/

// Alignment parameters that determine the cached comparisons (the k-mer cutoff is re-applied to
// the cached distances; qualities and error rates only enter lambda, which is recomputed).
struct CacheParams {
  int match, mismatch, gap, homo_gap, band_size, sse;
  bool use_kmers, vectorized_alignment, gapless;
  bool operator==(const CacheParams &o) const {
    return std::tie(match, mismatch, gap, homo_gap, band_size, sse, use_kmers, vectorized_alignment, gapless) ==
           std::tie(o.match, o.mismatch, o.gap, o.homo_gap, o.band_size, o.sse, o.use_kmers, o.vectorized_alignment, o.gapless);
  }
};

struct SubDelete {
  void operator()(Sub *s) const { sub_free(s); }
};

struct CompCacheSlot {
  dada2::AlignmentCache::Impl *owner = nullptr;
  bool valid = false;
  CacheParams params{};
  std::vector<std::string> seqs;  // the sample the slot belongs to (compared exactly)
  std::unordered_map<unsigned int, std::unique_ptr<CentreCache>> centres;
  // Final substitutions (centre, raw) -> Sub of earlier runs (computed without the k-mer screen,
  // so they can differ from the b_compare alignments and are cached separately).
  std::unordered_map<uint64_t, std::unique_ptr<Sub, SubDelete>> final_subs;
  size_t bytes = 0;               // budget used by this slot

  // Move-only, declared explicitly: MSVC's unordered_map move is not noexcept, so without the
  // deleted copy std::vector::resize would instantiate the (ill-formed) implicit copy.
  CompCacheSlot() = default;
  CompCacheSlot(const CompCacheSlot &) = delete;
  CompCacheSlot &operator=(const CompCacheSlot &) = delete;
  CompCacheSlot(CompCacheSlot &&) = default;
  CompCacheSlot &operator=(CompCacheSlot &&) = default;
};

struct dada2::AlignmentCache::Impl {
  std::vector<CompCacheSlot> slots;
  size_t budget = 0;
  std::atomic<size_t> used{0};
};

bool cache_reserve(CompCacheSlot *slot, size_t bytes) {
  auto &used = slot->owner->used;
  size_t cur = used.load();
  do {
    if(cur + bytes > slot->owner->budget) { return false; }
  } while(!used.compare_exchange_weak(cur, cur + bytes));
  slot->bytes += bytes;
  return true;
}

CentreCache *cache_centre(CompCacheSlot *slot, unsigned int centre, unsigned int nraw) {
  auto it = slot->centres.find(centre);
  if(it != slot->centres.end()) { return it->second.get(); }
  if(!cache_reserve(slot, (size_t)nraw * (sizeof(double) + sizeof(uint16_t) + sizeof(uint32_t)) + 128)) { return NULL; }
  auto cc = std::make_unique<CentreCache>();
  cc->kdist.assign(nraw, 0.0);
  cc->nsubs.assign(nraw, CACHE_UNKNOWN);
  cc->off.assign(nraw, 0);
  CentreCache *ptr = cc.get();
  slot->centres.emplace(centre, std::move(cc));
  return ptr;
}

namespace {

void slot_reset(CompCacheSlot &slot) {
  slot.owner->used.fetch_sub(slot.bytes);
  slot.bytes = 0;
  slot.centres.clear();
  slot.final_subs.clear();
  slot.seqs.clear();
  slot.valid = false;
}

// The slot for this sample and these parameters (reset if it held another sample), or NULL if
// the sample does not fit in the budget.
CompCacheSlot *slot_for(dada2::AlignmentCache *cache, size_t index, const std::vector<std::string> &seqs, const CacheParams &params) {
  if(!cache) { return NULL; }
  auto &impl = cache->impl();
  if(index >= impl.slots.size()) { dada2_stop_invalid("Alignment cache slot out of range."); }
  CompCacheSlot &slot = impl.slots[index];
  if(slot.valid && slot.params == params && slot.seqs == seqs) { return &slot; }
  slot_reset(slot);
  size_t bytes = 0;
  for(const auto &sq : seqs) { bytes += sq.size() + sizeof(std::string); }
  if(!cache_reserve(&slot, bytes)) { return NULL; }
  slot.seqs = seqs;
  slot.params = params;
  slot.valid = true;
  return &slot;
}

}  // namespace

namespace dada2 {

AlignmentCache::AlignmentCache(size_t n_slots, size_t budget_bytes) : impl_(std::make_unique<Impl>()) {
  impl_->budget = budget_bytes;
  impl_->slots.resize(n_slots);
  for(auto &slot : impl_->slots) { slot.owner = impl_.get(); }
}
AlignmentCache::~AlignmentCache() = default;
size_t AlignmentCache::n_slots() const { return impl_->slots.size(); }
size_t AlignmentCache::budget_bytes() const { return impl_->budget; }
size_t AlignmentCache::bytes_used() const { return impl_->used.load(); }

//------------------------------------------------------------------
// C interface to run DADA on the provided unique sequences/abundance pairs (was the Rcpp
// export dada_uniques). err: 16 x err_ncol row-major; quals: nraw x quals_ncol row-major.
DadaOutput dada_uniques(const std::vector<std::string> &seqs, const std::vector<int> &abundances,
                        const std::vector<bool> &priors, const double *err, int err_ncol,
                        const double *quals, int quals_ncol, const DadaOptions &opt, int n_threads,
                        bool verbose, AlignmentCache *cache, size_t cache_slot) {
  unsigned int i, r, index, pos, nraw, maxlen, minlen;
  Raw *raw;
  int SSE = opt.sse;

  /********** INPUT VALIDATION *********/
  // Check lengths of seqs and abundances vectors
  nraw = seqs.size();
  if(nraw == 0) { dada2_stop_invalid("Zero input sequences."); }
  if(nraw != abundances.size()) { dada2_stop_invalid("Sequence and abundance vectors had different lengths."); }
  if(nraw != priors.size()) { dada2_stop_invalid("Sequence and priors vectors had different lengths."); }
  maxlen=0;
  minlen=SEQLEN;
  for(index=0;index<nraw;index++) {
    if(seqs[index].length() > maxlen) { maxlen = seqs[index].length(); }
    if(seqs[index].length() < minlen) { minlen = seqs[index].length(); }
  }
  if(maxlen >= SEQLEN) { dada2_stop_invalid("Input sequences exceed the maximum allowed string length."); }
  if(minlen <= KMER_SIZE) { dada2_stop_invalid("Input sequences must all be longer than the kmer-size (%i).", KMER_SIZE); }
  for(index=0;index<nraw;index++) {
    if(abundances[index] < 1) { dada2_stop_invalid("Abundances must be positive integers."); }
  }
  for(bool ok : is_acgt(seqs)) {
    if(!ok) { dada2_stop_invalid("Sequences must contain only A/C/G/T characters."); }
  }

  // Check for presence of quality scores and their lengths
  bool has_quals = false;
  if(quals && quals_ncol > 0) { // Each sequence is a ROW, each column is a POSITION
    has_quals = true;
    if((unsigned int)quals_ncol != maxlen) {
      dada2_stop_invalid("Sequence must have associated qualities for each nucleotide position.");
    }
  }
  // Check error matrix
  if(err == nullptr || err_ncol < 1) {
    dada2_stop_invalid("Error matrix must have 16 rows and at least one column.");
  }
  for(size_t e=0;e<16*(size_t)err_ncol;e++) {
    if(!(err[e] >= 0.0 && err[e] <= 1.0)) { dada2_stop_invalid("Error matrix entries must be finite and in [0, 1]."); }
  }
  if(opt.use_quals && !has_quals) {
    dada2_stop_invalid("Quality scores are required when use_quals is true.");
  }
  ErrMat errMat{err, (unsigned int)err_ncol};
  if(has_quals) {
    for(index=0;index<nraw;index++) {
      for(pos=0;pos<seqs[index].length();pos++) {
        double q = quals[(size_t)index*quals_ncol + pos];
        if(!(q >= 0) || round(q) > 255) { dada2_stop_invalid("Invalid quality score (must be in 0..255, no NaN)."); }
      }
    }
  }
  // pydada2: upstream set SSE=0 on non-x86; kmers.cpp provides exact scalar emulations of the
  // SSE2 routines instead, so SSE modes give x86-identical results on every platform.

  DadaState st;
  /********** CONSTRUCT RAWS *********/
  st.raws.assign(nraw, nullptr);
  Raw **raws = st.raws.data();
  // Construct a raw for each input sequence, store in raws[index]
  dada2::parallel_for(0, nraw, [&](std::size_t begin, std::size_t end) {
    std::vector<char> seq;
    for(std::size_t index=begin;index<end;index++) {
      seq.assign(seqs[index].c_str(), seqs[index].c_str() + seqs[index].size() + 1);
      nt2int(seq.data(), seq.data());
      double *qual = has_quals ? const_cast<double *>(quals + index*(size_t)quals_ncol) : NULL;
      raws[index] = raw_new(seq.data(), qual, abundances[index], priors[index]);
      raws[index]->index = index;
    }
  }, 64, n_threads);

  if(opt.use_kmers) {
    size_t n_kmer = 1 << (2*KMER_SIZE);
    st.k8.assign(nraw * n_kmer, 0);
    st.k16.assign(nraw * n_kmer, 0);
    st.kord.assign((size_t)nraw * maxlen + 8, 0); // +8: SSE loads read whole 8-element blocks
    dada2::parallel_for(0, nraw, [&](std::size_t begin, std::size_t end) {
      for(std::size_t index=begin;index<end;index++) {
        Raw *raw = raws[index];
        /// Add uint8_t kmer index in contiguous memory block
        raw->kmer8 = &st.k8[index*n_kmer];
        assign_kmer8(raw->kmer8, raw->seq, KMER_SIZE);
        /// Add uint16_t kmer index in contiguous memory block
        raw->kmer = &st.k16[index*n_kmer];
        assign_kmer(raw->kmer, raw->seq, KMER_SIZE);
        /// Add uint16_t ordered kmer record in contiguous memory block
        raw->kord = &st.kord[index*(size_t)maxlen];
        assign_kmer_order(raw->kord, raw->seq, KMER_SIZE);
      }
    }, 64, n_threads);
  } else {
    for(index=0;index<nraw;index++) {
      raw = raws[index];
      raw->kmer8 = NULL;
      raw->kmer = NULL;
      raw->kord = NULL;
    }
  }

  int match = opt.match, mismatch = opt.mismatch, gap = opt.gap, homo_gap = opt.homo_gap;
  int band_size = opt.band_size;
  bool use_kmers = opt.use_kmers, vectorized_alignment = opt.vectorized_alignment, gapless = opt.gapless;

  /********** RUN DADA *********/
  st.bb = b_new(raws, nraw, opt.omega_a, opt.omega_p, opt.use_quals); // New cluster with all sequences in 1 bi
  B *bb = st.bb;
  bb->cache = slot_for(cache, cache_slot, seqs, CacheParams{match, mismatch, gap, homo_gap, band_size, SSE, use_kmers, vectorized_alignment, gapless});
  LambdaTable lambda_table;  // must outlive run_dada
  lambda_table_init(lambda_table, raws, nraw, maxlen, errMat, opt.use_quals, n_threads);
  bb->lam = &lambda_table;
  run_dada(bb, errMat,
           match, mismatch, gap, homo_gap, use_kmers, opt.kdist_cutoff, band_size,
           opt.detect_singletons, opt.max_clust, opt.min_fold, opt.min_hamming, opt.min_abund,
           vectorized_alignment, n_threads, verbose, SSE, gapless, opt.greedy);

  /********** MAKE OUTPUT *********/
  // Create subs for all the relevant alignments
  st.subs.assign(bb->nraw, nullptr);
  st.birth_subs.assign(bb->nclust, nullptr);
  Sub **subs = st.subs.data();
  Sub **birth_subs = st.birth_subs.data();

  // Was FinalSubsParallel (over clusters); now over raws, then birth subs over clusters. Each
  // raw/cluster writes only its own slot. With a comparison cache, the final substitutions of
  // earlier runs on the same sample are reused (same alignment call, same inputs).
  std::vector<std::pair<unsigned int, Raw *>> members;
  members.reserve(bb->nraw);
  for(i=0;i<bb->nclust;i++) {
    for(r=0;r<bb->bi[i]->nraw;r++) { members.emplace_back(i, bb->bi[i]->raw[r]); }
  }
  CompCacheSlot *slot = bb->cache;
  auto final_key = [&](std::size_t m) {
    return ((uint64_t)bb->bi[members[m].first]->center->index << 32) | members[m].second->index;
  };
  std::vector<unsigned char> from_cache(slot ? members.size() : 0, 0);
  dada2::parallel_for(0, members.size(), [&](std::size_t begin, std::size_t end) {
    for(std::size_t m=begin;m<end;m++) {
      Raw *raw = members[m].second;
      if(slot) { // the alignment of an earlier run (same sample and parameters)
        auto it = slot->final_subs.find(final_key(m));
        if(it != slot->final_subs.end()) {
          subs[raw->index] = sub_copy(it->second.get());
          from_cache[m] = 1;
          continue;
        }
      }
      subs[raw->index] = sub_new(bb->bi[members[m].first]->center, raw, match, mismatch, gap, homo_gap, false, 1.0, band_size, vectorized_alignment, SSE, gapless);
    }
  }, GRAIN_SIZE, n_threads);
  if(slot) {
    for(std::size_t m=0;m<members.size();m++) {
      Sub *sub = subs[members[m].second->index];
      if(from_cache[m] || !sub) { continue; }
      const size_t bytes = sizeof(Sub) + 64 + sub->len0 * sizeof(uint16_t) + sub->nsubs * (sizeof(uint16_t) + 4);
      if(!cache_reserve(slot, bytes)) { break; }
      slot->final_subs.emplace(final_key(m), std::unique_ptr<Sub, SubDelete>(sub_copy(sub)));
    }
  }
  dada2::parallel_for(0, bb->nclust, [&](std::size_t begin, std::size_t end) {
    for(std::size_t i=begin;i<end;i++) {
      // Make birth sub for that cluster
      if(i==0) { birth_subs[i] = NULL; }
      else {
        birth_subs[i] = sub_new(bb->bi[bb->bi[i]->birth_comp.i]->center, bb->bi[i]->center, match, mismatch, gap, homo_gap, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless);
      }
    }
  }, GRAIN_SIZE, n_threads);

  DadaOutput out;
  // Assign final w/in cluster P
  out.pval_raw.assign(nraw, 0.0);
  for(i=0;i<bb->nclust;i++) {
    for(r=0;r<bb->bi[i]->nraw;r++) {
      raw = bb->bi[i]->raw[r];
      index = raw->index;
      if(bb->bi[i]->center == raw) { // centers always have correct=true
        raw->p = 1.0; // Due to fp precision, possible to get 1's from non-centers in below calc
      } else {
        raw->p = calc_pA(raw->reads, raw->comp.lambda * bb->bi[i]->reads, true);
        if(raw->p < opt.omega_c) { raw->correct = false; } // Don't correct!
      }
      out.pval_raw[index] = raw->p;
    }
  }

  b_make_clustering_df(bb, subs, birth_subs, has_quals, out);
  b_make_transition_by_quality_matrix(bb, subs, has_quals, err_ncol, out);
  b_make_cluster_quality_matrix(bb, subs, has_quals, maxlen, out);
  b_make_birth_subs_df(bb, birth_subs, has_quals, out);

  // Make map from uniques to cluster
  out.map.assign(nraw, -1);
  for(i=0;i<bb->nclust;i++) {
    for(r=0;r<bb->bi[i]->nraw;r++) {
      raw = bb->bi[i]->raw[r];
      if(raw->correct) {
        out.map[raw->index] = i; // 0-based (upstream: i+1 for R)
      } else {
        out.map[raw->index] = -1; // upstream: NA_INTEGER
      }
    }
  }
  // DadaState releases subs, B and raws.
  return out;
}

//------------------------------------------------------------------
// Several samples, run concurrently on the shared pool (see api.h).
std::vector<DadaOutput> dada_uniques_batch(const std::vector<DadaInput> &samples,
                                           const DadaOptions &opt, int n_threads,
                                           int sample_threads, bool verbose,
                                           AlignmentCache *cache) {
  const size_t nsam = samples.size();
  std::vector<DadaOutput> out(nsam);
  if(nsam == 0) { return out; }
  if(cache && cache->n_slots() < nsam) { dada2_stop_invalid("Alignment cache has fewer slots than samples."); }
  for(const auto &in : samples) {
    if(!in.seqs || !in.abundances || !in.priors) { dada2_stop_invalid("Incomplete sample input."); }
  }
  // Largest samples first, so that the last ones to finish are small.
  std::vector<size_t> order(nsam);
  for(size_t k=0;k<nsam;k++) { order[k] = k; }
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return samples[a].seqs->size() > samples[b].seqs->size();
  });

  struct Batch {
    const std::vector<DadaInput> *samples;
    std::vector<DadaOutput> *out;
    const std::vector<size_t> *order;
    const DadaOptions *opt;
    AlignmentCache *cache;
    int n_threads;
    bool verbose;
    std::thread::id caller;
    std::atomic<bool> cancel{false};
    std::mutex m;
    std::condition_variable cv;
    size_t next = 0, running = 0;  // guarded by m
    std::exception_ptr error;      // guarded by m
    void fail(std::exception_ptr e) {
      std::lock_guard<std::mutex> lk(m);
      if(!error) { error = e; }
      cancel.store(true);
    }
  } batch;
  batch.samples = &samples;
  batch.out = &out;
  batch.order = &order;
  batch.opt = &opt;
  batch.cache = cache;
  batch.n_threads = n_threads;
  batch.verbose = verbose;
  batch.caller = std::this_thread::get_id();

  auto run = [](void *p) {
    Batch &bt = *static_cast<Batch *>(p);
    ScopedCancel cancel_scope(&bt.cancel);
    for(;;) {
      size_t s;
      {
        std::lock_guard<std::mutex> lk(bt.m);
        if(bt.cancel.load() || bt.next >= bt.order->size()) { break; }
        s = (*bt.order)[bt.next++];
        bt.running++;
      }
      try {
        const DadaInput &in = (*bt.samples)[s];
        (*bt.out)[s] = dada_uniques(*in.seqs, *in.abundances, *in.priors, in.err, in.err_ncol,
                                    in.quals, in.quals_ncol, *bt.opt, bt.n_threads, bt.verbose,
                                    bt.cache, s);
      } catch(...) {
        bt.fail(std::current_exception());
      }
      {
        std::lock_guard<std::mutex> lk(bt.m);
        bt.running--;
      }
      bt.cv.notify_all();
    }
    if(std::this_thread::get_id() != bt.caller) { return; }
    // The caller waits for the samples still running on other threads and keeps checking for
    // interrupts meanwhile (the check needs the caller's handler, e.g. Python signals).
    std::unique_lock<std::mutex> lk(bt.m);
    while(bt.running > 0) {
      if(bt.cv.wait_for(lk, std::chrono::milliseconds(50), [&] { return bt.running == 0; })) { break; }
      lk.unlock();
      try {
        dada2_check_interrupt();
      } catch(...) {
        bt.fail(std::current_exception());
      }
      lk.lock();
    }
  };
  const size_t workers = std::max<size_t>(1, std::min<size_t>(sample_threads < 1 ? 1 : (size_t)sample_threads, nsam));
  detail::run_on_pool(run, &batch, std::min(workers, detail::max_threads()) - 1);
  if(batch.error) { std::rethrow_exception(batch.error); }
  return out;
}

}  // namespace dada2

void run_dada(B *bb, const ErrMat &errMat,
              int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size,
              bool detect_singletons, int max_clust, double min_fold, int min_hamming, int min_abund,
              bool vectorized_alignment, int n_threads, bool verbose, int SSE, bool gapless, bool greedy) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;

  // Everyone gets aligned within the initial cluster, no KMER screen
  b_compare(bb, 0, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless, greedy, verbose, n_threads);
  b_p_update(bb, greedy, detect_singletons); // Calculates abundance p-value for each raw in its cluster (consensuses)
  
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, min_abund, verbose)) ) {
    if(verbose) dada2_log("\nNew Cluster C%i:", newi);
    b_compare(bb, newi, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless, greedy, verbose, n_threads);
    // Keep shuffling and updating until no more shuffles
    nshuffle = 0;
    do {
      shuffled = b_shuffle2(bb);
      if(verbose) { dada2_log("S"); }
    } while(shuffled && ++nshuffle < MAX_SHUFFLE);
    if(verbose && nshuffle >= MAX_SHUFFLE) { dada2_log("Warning: Reached maximum (%i) shuffles.\n", MAX_SHUFFLE); }

    b_p_update(bb, greedy, detect_singletons);
    dada2_check_interrupt();
  } // while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, min_abund, verbose)) )
  
  if(verbose) dada2_log("\nALIGN: %i aligns, %i shrouded (%i raw).\n", bb->nalign, bb->nshroud, bb->nraw);
}
