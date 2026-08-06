// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>
#include <vector>

#include "mapped_matrix.h"
#include "types.h"

namespace pgensparsescore {

struct BlockedDenseScorerStats {
  uint64_t block_ct = 0;
  uint32_t maximum_variant_ct = 0;
  uint64_t maximum_edge_ct = 0;
  uint64_t scoring_nanoseconds = 0;
};

// Buffers dense genotype rows, transposes their variant-major weight edges to
// score-major rows, and assigns each output score row to exactly one worker.
// This retains the output row in cache while all weights in the genotype block
// are applied and amortizes worker dispatch across many variants.
class BlockedDenseScorer {
 public:
  BlockedDenseScorer(uint32_t thread_ct, uint32_t score_ct,
                     uint32_t sample_ct, MappedMatrix* matrix);
  ~BlockedDenseScorer();

  BlockedDenseScorer(const BlockedDenseScorer&) = delete;
  BlockedDenseScorer& operator=(const BlockedDenseScorer&) = delete;

  void Add(const double* dosages, const std::vector<Edge>& edges);
  void Flush();
  const BlockedDenseScorerStats& stats() const { return stats_; }

 private:
  struct BlockEdge {
    uint32_t variant_idx;
    double beta;
  };

  void ApplyScore(uint32_t active_score_idx);
  void ApplyAvailableScores();
  void WorkerLoop();

  uint32_t thread_ct_;
  uint32_t sample_ct_;
  uint32_t variant_capacity_;
  MappedMatrix* matrix_;
  std::vector<double> dosages_;
  std::vector<std::vector<BlockEdge>> score_edges_;
  std::vector<uint32_t> active_scores_;
  uint64_t edge_ct_ = 0;

  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable start_;
  std::condition_variable done_;
  bool stopping_ = false;
  uint64_t generation_ = 0;
  uint32_t pending_worker_ct_ = 0;
  std::atomic<uint32_t> next_score_idx_{0};
  BlockedDenseScorerStats stats_;
};

}  // namespace pgensparsescore
