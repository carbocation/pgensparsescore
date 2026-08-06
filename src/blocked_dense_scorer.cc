// SPDX-License-Identifier: GPL-3.0-only
#include "blocked_dense_scorer.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace pgensparsescore {
namespace {

constexpr uint64_t kTargetDosageBufferBytes = 64ULL * 1024 * 1024;

void AddScaledDosages(const double* dosages, double beta, uint32_t sample_ct,
                      double* output) {
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    output[sample_idx] += beta * dosages[sample_idx];
  }
}

}  // namespace

BlockedDenseScorer::BlockedDenseScorer(uint32_t thread_ct, uint32_t score_ct,
                                       uint32_t sample_ct,
                                       MappedMatrix* matrix)
    : thread_ct_(thread_ct),
      sample_ct_(sample_ct),
      matrix_(matrix),
      score_edges_(score_ct) {
  if (!thread_ct_ || !score_ct || !sample_ct_ || !matrix_ ||
      matrix_->row_ct() != score_ct || matrix_->column_ct() != sample_ct_) {
    throw std::runtime_error("blocked dense scorer has incompatible shapes");
  }
  const uint64_t bytes_per_variant =
      static_cast<uint64_t>(sample_ct_) * sizeof(double);
  variant_capacity_ = static_cast<uint32_t>(std::max<uint64_t>(
      1, std::min<uint64_t>(UINT32_MAX,
                            kTargetDosageBufferBytes / bytes_per_variant)));
  dosages_.reserve(static_cast<uint64_t>(variant_capacity_) * sample_ct_);
  active_scores_.reserve(score_ct);
  workers_.reserve(thread_ct_ - 1);
  for (uint32_t worker_idx = 1; worker_idx < thread_ct_; ++worker_idx) {
    workers_.emplace_back(&BlockedDenseScorer::WorkerLoop, this);
  }
}

BlockedDenseScorer::~BlockedDenseScorer() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    stopping_ = true;
    ++generation_;
  }
  start_.notify_all();
  for (auto& worker : workers_) worker.join();
}

void BlockedDenseScorer::Add(const double* dosages,
                             const std::vector<Edge>& edges) {
  if (!dosages || edges.empty()) {
    throw std::runtime_error("cannot add an empty dense scoring variant");
  }
  const uint32_t variant_idx =
      static_cast<uint32_t>(dosages_.size() / sample_ct_);
  if (variant_idx == variant_capacity_) {
    Flush();
  }
  const uint32_t local_variant_idx =
      static_cast<uint32_t>(dosages_.size() / sample_ct_);
  const size_t old_size = dosages_.size();
  dosages_.resize(old_size + sample_ct_);
  std::memcpy(dosages_.data() + old_size, dosages,
              static_cast<size_t>(sample_ct_) * sizeof(double));
  for (const auto& edge : edges) {
    if (edge.score_idx >= score_edges_.size()) {
      throw std::runtime_error("dense block edge has invalid score index");
    }
    auto& score = score_edges_[edge.score_idx];
    if (score.empty()) active_scores_.push_back(edge.score_idx);
    score.push_back({local_variant_idx, edge.beta_alt});
  }
  edge_ct_ += edges.size();
}

void BlockedDenseScorer::ApplyScore(uint32_t active_score_idx) {
  const uint32_t score_idx = active_scores_[active_score_idx];
  double* output = matrix_->Row(score_idx);
  for (const auto& edge : score_edges_[score_idx]) {
    const double* input =
        dosages_.data() + static_cast<uint64_t>(edge.variant_idx) * sample_ct_;
    AddScaledDosages(input, edge.beta, sample_ct_, output);
  }
}

void BlockedDenseScorer::ApplyAvailableScores() {
  for (;;) {
    const uint32_t active_score_idx = next_score_idx_.fetch_add(1);
    if (active_score_idx >= active_scores_.size()) return;
    ApplyScore(active_score_idx);
  }
}

void BlockedDenseScorer::WorkerLoop() {
  uint64_t observed_generation = 0;
  for (;;) {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      start_.wait(lock, [&] {
        return stopping_ || generation_ != observed_generation;
      });
      if (stopping_) return;
      observed_generation = generation_;
    }
    ApplyAvailableScores();
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!--pending_worker_ct_) done_.notify_one();
    }
  }
}

void BlockedDenseScorer::Flush() {
  if (dosages_.empty()) return;
  const uint32_t variant_ct =
      static_cast<uint32_t>(dosages_.size() / sample_ct_);
  ++stats_.block_ct;
  stats_.maximum_variant_ct =
      std::max(stats_.maximum_variant_ct, variant_ct);
  stats_.maximum_edge_ct = std::max(stats_.maximum_edge_ct, edge_ct_);
  const auto scoring_started = std::chrono::steady_clock::now();

  if (thread_ct_ == 1 || active_scores_.size() == 1) {
    for (uint32_t active_idx = 0; active_idx < active_scores_.size();
         ++active_idx) {
      ApplyScore(active_idx);
    }
  } else {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      next_score_idx_.store(0);
      pending_worker_ct_ = thread_ct_ - 1;
      ++generation_;
    }
    start_.notify_all();
    ApplyAvailableScores();
    std::unique_lock<std::mutex> lock(mutex_);
    done_.wait(lock, [&] { return pending_worker_ct_ == 0; });
  }
  stats_.scoring_nanoseconds += static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - scoring_started)
          .count());

  for (const uint32_t score_idx : active_scores_) {
    score_edges_[score_idx].clear();
  }
  active_scores_.clear();
  dosages_.clear();
  edge_ct_ = 0;
}

}  // namespace pgensparsescore
