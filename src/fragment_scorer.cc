// SPDX-License-Identifier: GPL-3.0-only
#include "fragment_scorer.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <filesystem>
#include <exception>
#include <iostream>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include "io.h"
#include "scorer.h"

namespace pgensparsescore {
namespace {

using Header = std::unordered_map<std::string, size_t>;

Header MakeHeader(const std::vector<std::string>& fields,
                  const std::string& path) {
  Header result;
  for (size_t idx = 0; idx < fields.size(); ++idx) {
    std::string name = fields[idx];
    if (!name.empty() && name.front() == '#') name.erase(0, 1);
    if (!result.emplace(name, idx).second) {
      throw std::runtime_error(path + " has duplicate column " + name);
    }
  }
  return result;
}

size_t RequireColumn(const Header& header, const std::string& name,
                     const std::string& path) {
  const auto iter = header.find(name);
  if (iter == header.end()) {
    throw std::runtime_error(path + " is missing column " + name);
  }
  return iter->second;
}

enum class TileVariantStorage : uint8_t {
  kMissingVariant,
  kOmittedFrequency,
  kDense,
  kSparse,
};

struct TileVariantState {
  TileVariantStorage storage = TileVariantStorage::kMissingVariant;
  uint32_t storage_idx = 0;
};

struct StoredSparseDosage {
  double common = 0.0;
  double mean = 0.0;
  std::vector<uint32_t> sample_ids;
  std::vector<uint16_t> dosage16;
};

struct ScoreRowTask {
  const ScoreFragmentScoreRow* row = nullptr;
  uint32_t output_score_idx = 0;
};

struct ScoreRowWorkStats {
  uint64_t row_ct = 0;
  uint64_t scored_edge_ct = 0;
  uint64_t omitted_frequency_edge_ct = 0;
  uint64_t dense_edge_ct = 0;
  uint64_t sparse_edge_ct = 0;
  uint64_t dense_update_ct = 0;
  uint64_t sparse_update_ct = 0;
};

struct DenseEdge {
  uint32_t output_score_idx = 0;
  uint32_t dense_variant_idx = 0;
  double beta_alt = 0.0;
};

struct BlockedDenseWeight {
  uint32_t dense_variant_idx = 0;
  double beta_alt = 0.0;
};

struct BlockedDenseRow {
  double* output = nullptr;
  uint32_t weight_begin = 0;
  uint32_t weight_end = 0;
};

struct BlockedDensePlan {
  std::vector<BlockedDenseRow> rows;
  std::vector<BlockedDenseWeight> weights;
};

constexpr uint32_t kDenseSampleBlockSize = 256;

void AddScaledDosages(const double* __restrict__ input, double beta,
                      uint32_t sample_ct, double* __restrict__ output) {
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    output[sample_idx] += beta * input[sample_idx];
  }
}

class ScoreRowWorkers {
 public:
  ScoreRowWorkers(uint32_t thread_ct, uint32_t sample_ct,
                  std::vector<double>* baselines, Catalog* catalog,
                  MappedMatrix* matrix)
      : thread_ct_(thread_ct),
        sample_ct_(sample_ct),
        baselines_(baselines),
        catalog_(catalog),
        matrix_(matrix),
        worker_stats_(thread_ct),
        dense_edges_(thread_ct) {
    if (!thread_ct_ || !sample_ct_ || !baselines_ || !catalog_ || !matrix_) {
      throw std::runtime_error("invalid score-major worker configuration");
    }
    workers_.reserve(thread_ct_ - 1);
    for (uint32_t worker_idx = 1; worker_idx < thread_ct_; ++worker_idx) {
      workers_.emplace_back(&ScoreRowWorkers::WorkerLoop, this, worker_idx);
    }
  }

  ~ScoreRowWorkers() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
      ++generation_;
    }
    start_.notify_all();
    for (auto& worker : workers_) worker.join();
  }

  ScoreRowWorkStats Dispatch(
      const std::vector<ScoreRowTask>& tasks,
      const std::vector<TileVariantState>& variants,
      const std::vector<StoredSparseDosage>& sparse_dosages) {
    if (tasks.empty()) return {};
    {
      std::lock_guard<std::mutex> lock(mutex_);
      tasks_ = &tasks;
      variants_ = &variants;
      sparse_dosages_ = &sparse_dosages;
      next_task_.store(0);
      worker_error_ = nullptr;
      std::fill(worker_stats_.begin(), worker_stats_.end(),
                ScoreRowWorkStats{});
      for (auto& edges : dense_edges_) edges.clear();
      pending_worker_ct_ = thread_ct_ - 1;
      ++generation_;
    }
    start_.notify_all();
    ApplyAvailable(0);
    if (thread_ct_ > 1) {
      std::unique_lock<std::mutex> lock(mutex_);
      done_.wait(lock, [&] { return pending_worker_ct_ == 0; });
    }
    if (worker_error_) std::rethrow_exception(worker_error_);
    ScoreRowWorkStats result;
    for (const auto& input : worker_stats_) {
      result.row_ct += input.row_ct;
      result.scored_edge_ct += input.scored_edge_ct;
      result.omitted_frequency_edge_ct += input.omitted_frequency_edge_ct;
      result.dense_edge_ct += input.dense_edge_ct;
      result.sparse_edge_ct += input.sparse_edge_ct;
      result.dense_update_ct += input.dense_update_ct;
      result.sparse_update_ct += input.sparse_update_ct;
    }
    return result;
  }

  const std::vector<std::vector<DenseEdge>>& dense_edges() const {
    return dense_edges_;
  }

 private:
  void ApplyTask(const ScoreRowTask& task, uint32_t worker_idx,
                 ScoreRowWorkStats* stats) {
    if (!task.row || task.output_score_idx >= catalog_->scores.size()) {
      throw std::runtime_error("invalid score-major row task");
    }
    auto& score = catalog_->scores[task.output_score_idx];
    double* output = matrix_->Row(task.output_score_idx);
    uint32_t previous_variant_idx = 0;
    bool has_previous = false;
    for (uint32_t edge_idx = 0; edge_idx < task.row->edge_ct(); ++edge_idx) {
      const auto edge = task.row->edge(edge_idx);
      if (has_previous && edge.local_variant_idx < previous_variant_idx) {
        throw std::runtime_error("score-major row variants are out of order");
      }
      previous_variant_idx = edge.local_variant_idx;
      has_previous = true;
      const auto& variant = variants_->at(edge.local_variant_idx);
      if (variant.storage == TileVariantStorage::kMissingVariant) {
        ++score.missing_variant_ct;
        continue;
      }
      ++score.matched_weight_ct;
      if (edge.ref_effect) {
        ++score.ref_effect_ct;
      } else {
        ++score.alt_effect_ct;
      }
      if (variant.storage == TileVariantStorage::kOmittedFrequency) {
        ++score.missing_frequency_ct;
        ++stats->omitted_frequency_edge_ct;
        continue;
      }
      if (edge.ref_effect) {
        const double intercept = -2.0 * edge.beta_alt;
        (*baselines_)[task.output_score_idx] += intercept;
        catalog_->intercepts[task.output_score_idx] += intercept;
        score.ref_effect_intercept += intercept;
      }
      ++stats->scored_edge_ct;
      if (variant.storage == TileVariantStorage::kDense) {
        dense_edges_[worker_idx].push_back(
            {task.output_score_idx, variant.storage_idx, edge.beta_alt});
        ++stats->dense_edge_ct;
        stats->dense_update_ct += sample_ct_;
      } else {
        const auto& sparse = sparse_dosages_->at(variant.storage_idx);
        (*baselines_)[task.output_score_idx] += edge.beta_alt * sparse.common;
        for (uint32_t value_idx = 0; value_idx < sparse.sample_ids.size();
             ++value_idx) {
          const double dosage = sparse.dosage16[value_idx] == UINT16_MAX
                                    ? sparse.mean
                                    : static_cast<double>(
                                          sparse.dosage16[value_idx]) /
                                          16384.0;
          output[sparse.sample_ids[value_idx]] +=
              edge.beta_alt * (dosage - sparse.common);
        }
        ++stats->sparse_edge_ct;
        stats->sparse_update_ct += sparse.sample_ids.size();
      }
    }
    ++stats->row_ct;
  }

  void ApplyAvailable(uint32_t worker_idx) {
    try {
      for (;;) {
        const uint32_t task_idx = next_task_.fetch_add(1);
        if (task_idx >= tasks_->size()) return;
        ApplyTask((*tasks_)[task_idx], worker_idx, &worker_stats_[worker_idx]);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!worker_error_) worker_error_ = std::current_exception();
      next_task_.store(UINT32_MAX);
    }
  }

  void WorkerLoop(uint32_t worker_idx) {
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
      ApplyAvailable(worker_idx);
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!--pending_worker_ct_) done_.notify_one();
      }
    }
  }

  uint32_t thread_ct_;
  uint32_t sample_ct_;
  std::vector<double>* baselines_;
  Catalog* catalog_;
  MappedMatrix* matrix_;
  std::vector<std::thread> workers_;
  std::vector<ScoreRowWorkStats> worker_stats_;
  std::vector<std::vector<DenseEdge>> dense_edges_;
  std::mutex mutex_;
  std::condition_variable start_;
  std::condition_variable done_;
  bool stopping_ = false;
  uint64_t generation_ = 0;
  uint32_t pending_worker_ct_ = 0;
  std::atomic<uint32_t> next_task_{0};
  const std::vector<ScoreRowTask>* tasks_ = nullptr;
  const std::vector<TileVariantState>* variants_ = nullptr;
  const std::vector<StoredSparseDosage>* sparse_dosages_ = nullptr;
  std::exception_ptr worker_error_;
};

BlockedDensePlan BuildBlockedDensePlan(
    const std::vector<std::vector<DenseEdge>>& worker_edges,
    uint32_t dense_variant_ct, Catalog* catalog, MappedMatrix* matrix) {
  if (!catalog || !matrix || matrix->row_ct() != catalog->scores.size()) {
    throw std::runtime_error("invalid blocked-dense plan inputs");
  }
  std::vector<uint32_t> counts(catalog->scores.size(), 0);
  uint64_t weight_ct = 0;
  for (const auto& edges : worker_edges) {
    for (const auto& edge : edges) {
      if (edge.output_score_idx >= counts.size() ||
          edge.dense_variant_idx >= dense_variant_ct ||
          counts[edge.output_score_idx] == UINT32_MAX) {
        throw std::runtime_error("invalid blocked-dense edge");
      }
      ++counts[edge.output_score_idx];
      ++weight_ct;
    }
  }
  if (weight_ct > UINT32_MAX) {
    throw std::runtime_error("blocked-dense tile contains too many weights");
  }
  std::vector<uint32_t> offsets(counts.size() + 1, 0);
  for (uint32_t score_idx = 0; score_idx < counts.size(); ++score_idx) {
    offsets[score_idx + 1] = offsets[score_idx] + counts[score_idx];
  }
  BlockedDensePlan plan;
  plan.weights.resize(static_cast<size_t>(weight_ct));
  std::vector<uint32_t> next = offsets;
  for (const auto& edges : worker_edges) {
    for (const auto& edge : edges) {
      plan.weights[next[edge.output_score_idx]++] =
          {edge.dense_variant_idx, edge.beta_alt};
    }
  }
  plan.rows.reserve(catalog->scores.size());
  for (uint32_t score_idx = 0; score_idx < counts.size(); ++score_idx) {
    if (!counts[score_idx]) continue;
    plan.rows.push_back(
        {matrix->Row(score_idx), offsets[score_idx], offsets[score_idx + 1]});
  }
  return plan;
}

class BlockedDenseWorkers {
 public:
  BlockedDenseWorkers(uint32_t thread_ct, uint32_t sample_ct)
      : thread_ct_(thread_ct), sample_ct_(sample_ct) {
    if (!thread_ct_ || !sample_ct_) {
      throw std::runtime_error("invalid blocked-dense worker configuration");
    }
    workers_.reserve(thread_ct_ - 1);
    for (uint32_t worker_idx = 1; worker_idx < thread_ct_; ++worker_idx) {
      workers_.emplace_back(&BlockedDenseWorkers::WorkerLoop, this);
    }
  }

  ~BlockedDenseWorkers() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
      ++generation_;
    }
    start_.notify_all();
    for (auto& worker : workers_) worker.join();
  }

  uint32_t Dispatch(const BlockedDensePlan& plan,
                    const std::vector<double>& dense_dosages,
                    uint32_t dense_variant_ct) {
    if (plan.weights.empty()) return 0;
    if (!dense_variant_ct ||
        dense_dosages.size() !=
            static_cast<uint64_t>(dense_variant_ct) * sample_ct_) {
      throw std::runtime_error("blocked-dense dosage shape changed");
    }
    const uint32_t block_ct =
        (sample_ct_ + kDenseSampleBlockSize - 1) / kDenseSampleBlockSize;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      plan_ = &plan;
      dense_dosages_ = &dense_dosages;
      block_ct_ = block_ct;
      next_block_.store(0);
      worker_error_ = nullptr;
      pending_worker_ct_ = thread_ct_ - 1;
      ++generation_;
    }
    start_.notify_all();
    ApplyAvailable();
    if (thread_ct_ > 1) {
      std::unique_lock<std::mutex> lock(mutex_);
      done_.wait(lock, [&] { return pending_worker_ct_ == 0; });
    }
    if (worker_error_) std::rethrow_exception(worker_error_);
    return block_ct;
  }

 private:
  void ApplyBlock(uint32_t block_idx) const {
    const uint32_t sample_begin = block_idx * kDenseSampleBlockSize;
    const uint32_t sample_end =
        std::min(sample_ct_, sample_begin + kDenseSampleBlockSize);
    const uint32_t block_sample_ct = sample_end - sample_begin;
    for (const auto& row : plan_->rows) {
      double* output = row.output + sample_begin;
      for (uint32_t weight_idx = row.weight_begin;
           weight_idx < row.weight_end; ++weight_idx) {
        const auto& weight = plan_->weights[weight_idx];
        const double* input =
            dense_dosages_->data() +
            static_cast<uint64_t>(weight.dense_variant_idx) * sample_ct_ +
            sample_begin;
        AddScaledDosages(input, weight.beta_alt, block_sample_ct, output);
      }
    }
  }

  void ApplyAvailable() {
    try {
      for (;;) {
        const uint32_t block_idx = next_block_.fetch_add(1);
        if (block_idx >= block_ct_) return;
        ApplyBlock(block_idx);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!worker_error_) worker_error_ = std::current_exception();
      next_block_.store(UINT32_MAX);
    }
  }

  void WorkerLoop() {
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
      ApplyAvailable();
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!--pending_worker_ct_) done_.notify_one();
      }
    }
  }

  uint32_t thread_ct_ = 0;
  uint32_t sample_ct_ = 0;
  uint32_t block_ct_ = 0;
  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable start_;
  std::condition_variable done_;
  bool stopping_ = false;
  uint64_t generation_ = 0;
  uint32_t pending_worker_ct_ = 0;
  std::atomic<uint32_t> next_block_{0};
  const BlockedDensePlan* plan_ = nullptr;
  const std::vector<double>* dense_dosages_ = nullptr;
  std::exception_ptr worker_error_;
};

}  // namespace

std::vector<std::string> ReadScoreFragmentList(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) throw std::runtime_error(path + " is empty");
  const Header header = MakeHeader(SplitTabs(line), path);
  size_t fragment_idx = static_cast<size_t>(-1);
  for (const char* name : {"FRAGMENT", "PATH"}) {
    const auto iter = header.find(name);
    if (iter != header.end()) {
      fragment_idx = iter->second;
      break;
    }
  }
  if (fragment_idx == static_cast<size_t>(-1)) {
    throw std::runtime_error(path + " is missing FRAGMENT or PATH column");
  }
  const auto base = std::filesystem::absolute(path).parent_path();
  std::vector<std::string> result;
  std::unordered_set<std::string> seen;
  uint64_t line_number = 1;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty()) continue;
    const auto fields = SplitTabs(line);
    if (fields.size() <= fragment_idx || fields[fragment_idx].empty()) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has no fragment path");
    }
    std::filesystem::path resolved(fields[fragment_idx]);
    if (resolved.is_relative()) resolved = base / resolved;
    const std::string normalized = resolved.lexically_normal().string();
    if (!seen.insert(normalized).second) {
      throw std::runtime_error(path + " contains duplicate fragment " +
                               fields[fragment_idx]);
    }
    result.push_back(normalized);
  }
  if (result.empty()) throw std::runtime_error(path + " has no fragments");
  return result;
}

IndexedPvarStats AddIndexedPvar(
    const std::string& path, uint32_t input_idx, const VariantIndex& index,
    std::vector<IndexedVariantLocation>* locations,
    ProgressReporter* progress) {
  if (locations->size() != index.variant_ct()) {
    throw std::runtime_error("variant-location vector does not match index");
  }
  LineReader reader(path);
  std::string line;
  Header header;
  size_t id_idx = 0;
  size_t ref_idx = 0;
  size_t alt_idx = 0;
  size_t maximum = 0;
  uint64_t line_number = 0;
  IndexedPvarStats stats;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty() || line.rfind("##", 0) == 0) continue;
    if (header.empty()) {
      header = MakeHeader(SplitTabs(line), path);
      id_idx = RequireColumn(header, "ID", path);
      ref_idx = RequireColumn(header, "REF", path);
      alt_idx = RequireColumn(header, "ALT", path);
      maximum = std::max({id_idx, ref_idx, alt_idx});
      continue;
    }
    if (stats.row_ct == UINT32_MAX) {
      throw std::runtime_error(path + " exceeds the supported PVAR row count");
    }
    const uint32_t pgen_variant_idx = stats.row_ct++;
    const auto fields = SplitTabs(line);
    if (fields.size() <= maximum) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has too few fields");
    }
    const auto ordinal = index.Lookup(fields[id_idx]);
    if (!ordinal) {
      if (progress && !(stats.row_ct % 1000000)) {
        progress->MaybeEvent("score", "read_pvar",
                             {{"pgen_input_index", input_idx},
                              {"pvar_rows_scanned", stats.row_ct},
                              {"indexed_variants_matched", stats.matched_variant_ct}},
                             {{"pvar", path}});
      }
      continue;
    }
    if (fields[alt_idx].find(',') != std::string::npos ||
        fields[ref_idx] != index.ref(*ordinal) ||
        fields[alt_idx] != index.alt(*ordinal)) {
      throw std::runtime_error(
          path + ": line " + std::to_string(line_number) +
          " alleles disagree with the variant index for " + fields[id_idx]);
    }
    auto& location = locations->at(*ordinal);
    if (location.present()) {
      throw std::runtime_error(
          "indexed variant occurs in more than one PVAR row: " +
          fields[id_idx]);
    }
    location.input_idx = input_idx;
    location.pgen_variant_idx = pgen_variant_idx;
    ++stats.matched_variant_ct;
    if (progress && !(stats.row_ct % 1000000)) {
      progress->MaybeEvent("score", "read_pvar",
                           {{"pgen_input_index", input_idx},
                            {"pvar_rows_scanned", stats.row_ct},
                            {"indexed_variants_matched", stats.matched_variant_ct}},
                           {{"pvar", path}});
    }
  }
  if (header.empty()) throw std::runtime_error(path + " has no header");
  if (!stats.row_ct) throw std::runtime_error(path + " has no variants");
  if (progress) {
    progress->Event("score", "pvar_loaded",
                    {{"pgen_input_index", input_idx},
                     {"pvar_rows", stats.row_ct},
                     {"indexed_variants_matched", stats.matched_variant_ct}},
                    {{"pvar", path}});
  }
  return stats;
}

LoadedScoreFragments LoadScoreFragments(
    const std::vector<std::string>& paths, const VariantIndex& index,
    const std::string& score_schema_path, ProgressReporter* progress) {
  if (paths.empty()) throw std::runtime_error("no score fragments supplied");
  LoadedScoreFragments result;
  result.fragments.reserve(paths.size());
  uint64_t score_ct = 0;
  uint32_t tile_size = 0;
  uint32_t tile_ct = 0;
  for (const auto& path : paths) {
    result.fragments.emplace_back(path);
    const auto& fragment = result.fragments.back();
    if (fragment.variant_ct() != index.variant_ct() ||
        fragment.signature_lo() != index.signature_lo() ||
        fragment.signature_hi() != index.signature_hi()) {
      throw std::runtime_error(path +
                               " was built for a different variant index");
    }
    if (!tile_size) {
      tile_size = fragment.tile_size();
      tile_ct = fragment.tile_ct();
    } else if (fragment.tile_size() != tile_size ||
               fragment.tile_ct() != tile_ct) {
      throw std::runtime_error(
          "score fragments use different scoring-tile geometries");
    }
    score_ct += fragment.scores().size();
    result.weight_ct += fragment.weight_ct();
    result.file_byte_ct += fragment.file_bytes();
    if (progress) {
      progress->Event("score", "fragment_loaded",
                      {{"fragments_loaded", result.fragments.size()},
                       {"fragments_total", paths.size()},
                       {"scores_loaded", score_ct},
                       {"weights_available", result.weight_ct},
                       {"fragment_bytes", result.file_byte_ct}},
                      {{"fragment", path}});
    }
  }
  if (!score_ct || score_ct > UINT32_MAX) {
    throw std::runtime_error("score fragments contain too many scores");
  }

  std::vector<std::string> schema_ids;
  std::vector<std::string> schema_columns;
  if (!score_schema_path.empty()) {
    LineReader reader(score_schema_path);
    std::string line;
    if (!reader.GetLine(&line)) {
      throw std::runtime_error(score_schema_path + " is empty");
    }
    const Header header = MakeHeader(SplitTabs(line), score_schema_path);
    size_t id_idx = static_cast<size_t>(-1);
    for (const char* name : {"SCORE_ID", "SCORE", "COLUMN_NAME"}) {
      const auto iter = header.find(name);
      if (iter != header.end()) {
        id_idx = iter->second;
        break;
      }
    }
    if (id_idx == static_cast<size_t>(-1)) {
      throw std::runtime_error(
          score_schema_path + " is missing SCORE_ID, SCORE, or COLUMN_NAME");
    }
    const auto column_iter = header.find("COLUMN_NAME");
    const size_t column_idx =
        column_iter == header.end() ? id_idx : column_iter->second;
    const size_t maximum = std::max(id_idx, column_idx);
    std::unordered_set<std::string> seen_schema_ids;
    std::unordered_set<std::string> seen_schema_columns;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) continue;
      const auto fields = SplitTabs(line);
      if (fields.size() <= maximum || fields[id_idx].empty() ||
          fields[column_idx].empty() ||
          !seen_schema_ids.insert(fields[id_idx]).second ||
          !seen_schema_columns.insert(fields[column_idx]).second) {
        throw std::runtime_error(
            score_schema_path + ": line " + std::to_string(line_number) +
            " has an empty or duplicate score ID/column");
      }
      schema_ids.push_back(fields[id_idx]);
      schema_columns.push_back(fields[column_idx]);
    }
    if (schema_ids.size() != score_ct) {
      throw std::runtime_error(
          "score schema and fragments contain different score counts");
    }
  }

  std::unordered_map<std::string, uint32_t> output_by_id;
  if (!schema_ids.empty()) {
    output_by_id.reserve(schema_ids.size());
    for (uint32_t idx = 0; idx < schema_ids.size(); ++idx) {
      output_by_id.emplace(schema_ids[idx], idx);
    }
  }
  result.catalog.scores.resize(static_cast<size_t>(score_ct));
  result.catalog.intercepts.assign(static_cast<size_t>(score_ct), 0.0);
  result.score_maps.resize(result.fragments.size());
  std::vector<bool> output_seen(score_ct, false);
  std::unordered_set<std::string> fragment_ids;
  std::unordered_set<std::string> output_columns;
  uint32_t next_output_idx = 0;
  for (uint32_t fragment_idx = 0; fragment_idx < result.fragments.size();
       ++fragment_idx) {
    const auto& fragment = result.fragments[fragment_idx];
    auto& score_map = result.score_maps[fragment_idx];
    score_map.reserve(fragment.scores().size());
    for (const auto& fragment_score : fragment.scores()) {
      if (!fragment_ids.insert(fragment_score.score_id).second) {
        throw std::runtime_error(
            "score fragments contain duplicate stable score ID " +
            fragment_score.score_id);
      }
      uint32_t output_idx = next_output_idx++;
      std::string output_column = fragment_score.info.id;
      if (!schema_ids.empty()) {
        const auto schema = output_by_id.find(fragment_score.score_id);
        if (schema == output_by_id.end()) {
          throw std::runtime_error("score schema has no row for fragment score " +
                                   fragment_score.score_id);
        }
        output_idx = schema->second;
        output_column = schema_columns[output_idx];
      }
      if (output_seen[output_idx] ||
          !output_columns.insert(output_column).second) {
        throw std::runtime_error(
            "score fragments resolve to a duplicate output column");
      }
      output_seen[output_idx] = true;
      score_map.push_back(output_idx);
      result.catalog.scores[output_idx] = fragment_score.info;
      auto& score = result.catalog.scores[output_idx];
      score.id = output_column;
      score.ref_effect_intercept = 0.0;
    }
  }
  if (std::find(output_seen.begin(), output_seen.end(), false) !=
      output_seen.end()) {
    throw std::runtime_error("score schema contains a score absent from fragments");
  }
  if (progress) {
    progress->Event("score", "fragments_ready",
                    {{"fragments", result.fragments.size()},
                     {"scores", result.catalog.scores.size()},
                     {"weights_available", result.weight_ct},
                     {"fragment_bytes", result.file_byte_ct}});
  }
  return result;
}

ScoreRunStats ScoreFragments(
    const VariantIndex& index,
    const std::vector<ScoreFragmentReader>& fragments,
    const std::vector<std::vector<uint32_t>>& score_maps,
    const std::vector<IndexedVariantLocation>& locations,
    const IndexedFrequencyTable* frequencies,
    MissingFrequencyPolicy missing_frequency_policy,
    const std::vector<PgenDosageReader*>& readers, Catalog* catalog,
    MappedMatrix* matrix, uint32_t thread_ct, ProgressReporter* progress) {
  if (locations.size() != index.variant_ct() || fragments.empty() ||
      score_maps.size() != fragments.size() ||
      readers.empty() || matrix->row_ct() != catalog->scores.size() ||
      matrix->column_ct() != readers.front()->sample_ct()) {
    throw std::runtime_error("fragment scoring inputs have incompatible shapes");
  }
  for (const auto* reader : readers) {
    if (!reader || reader->sample_ct() != matrix->column_ct()) {
      throw std::runtime_error("fragment PGEN readers have different samples");
    }
  }
  if (frequencies && frequencies->alt_dosage_means.size() != index.variant_ct()) {
    throw std::runtime_error("indexed frequencies do not match variant index");
  }
  if (!frequencies && missing_frequency_policy != MissingFrequencyPolicy::kCohort) {
    throw std::runtime_error(
        "error and omit missing-frequency policies require frequencies");
  }
  const uint32_t tile_size = fragments.front().tile_size();
  const uint32_t tile_ct = fragments.front().tile_ct();
  for (const auto& fragment : fragments) {
    if (fragment.tile_size() != tile_size || fragment.tile_ct() != tile_ct) {
      throw std::runtime_error("score fragments use different scoring tiles");
    }
  }

  ScoreRunStats stats;
  for (const auto& score : catalog->scores) {
    stats.omitted_frequency_edge_ct += score.missing_frequency_ct;
  }
  std::vector<double> baselines(catalog->scores.size(), 0.0);
  ScoreRowWorkers scoring_workers(thread_ct, matrix->column_ct(), &baselines,
                                  catalog, matrix);
  BlockedDenseWorkers dense_workers(thread_ct, matrix->column_ct());
  stats.blocked_dense_sample_block_size = kDenseSampleBlockSize;
  uint64_t variant_groups_processed = 0;
  constexpr uint64_t kMaximumCopiedSparseBytes = 512ULL * 1024 * 1024;

  for (uint32_t tile_idx = 0; tile_idx < tile_ct; ++tile_idx) {
    std::vector<ScoreFragmentTile> tiles;
    tiles.reserve(fragments.size());
    uint32_t tile_variant_ct = 0;
    for (uint32_t fragment_idx = 0; fragment_idx < fragments.size();
         ++fragment_idx) {
      tiles.push_back(fragments[fragment_idx].OpenTile(tile_idx));
      if (!tile_variant_ct) {
        tile_variant_ct = tiles.back().variant_ct();
      } else if (tiles.back().variant_ct() != tile_variant_ct) {
        throw std::runtime_error("score-fragment tile sizes disagree");
      }
    }
    const uint32_t bitmap_word_ct = (tile_variant_ct + 63) / 64;
    std::vector<uint64_t> referenced_variants(bitmap_word_ct, 0);
    for (const auto& tile : tiles) {
      tile.OrReferencedVariants(&referenced_variants);
    }
    std::vector<TileVariantState> variant_states(tile_variant_ct);
    std::vector<double> dense_dosages;
    dense_dosages.reserve(static_cast<uint64_t>(tile_variant_ct) *
                          matrix->column_ct());
    std::vector<StoredSparseDosage> sparse_dosages;
    sparse_dosages.reserve(tile_variant_ct);
    uint64_t copied_sparse_bytes = 0;
    uint64_t referenced_in_tile = 0;
    const uint32_t first_ordinal = tile_idx * tile_size;
    for (uint32_t local_variant_idx = 0;
         local_variant_idx < tile_variant_ct; ++local_variant_idx) {
      if (!(referenced_variants[local_variant_idx / 64] &
            (uint64_t{1} << (local_variant_idx % 64)))) {
        continue;
      }
      ++referenced_in_tile;
      ++variant_groups_processed;
      const uint32_t ordinal = first_ordinal + local_variant_idx;
      const auto& location = locations[ordinal];
      if (!location.present()) {
        continue;
      }
      std::optional<double> imputation_mean;
      bool omit = false;
      if (frequencies) {
        const double mean = frequencies->alt_dosage_means[ordinal];
        if (std::isnan(mean)) {
          ++stats.missing_frequency_variant_ct;
          if (missing_frequency_policy == MissingFrequencyPolicy::kError) {
            throw std::runtime_error(
                "frequency file has no usable row for indexed variant ordinal " +
                std::to_string(ordinal));
          }
          if (missing_frequency_policy == MissingFrequencyPolicy::kOmit) {
            omit = true;
            ++stats.omitted_frequency_variant_ct;
          } else {
            ++stats.cohort_frequency_variant_ct;
          }
        } else {
          imputation_mean = mean;
          ++stats.external_frequency_variant_ct;
        }
      } else {
        ++stats.cohort_frequency_variant_ct;
      }
      if (omit) {
        variant_states[local_variant_idx].storage =
            TileVariantStorage::kOmittedFrequency;
        continue;
      }
      if (location.input_idx >= readers.size()) {
        throw std::runtime_error("indexed PVAR location has invalid input index");
      }
      PgenDosageReader* reader = readers[location.input_idx];
      const DosageView dosage =
          reader->Read(location.pgen_variant_idx, imputation_mean);
      ++stats.variant_ct;
      stats.imputed_value_ct += dosage.missing_ct;
      if (dosage.sparse) {
        ++stats.sparse_variant_ct;
        stats.sparse_value_ct += dosage.sparse_value_ct;
        const uint64_t value_bytes =
            static_cast<uint64_t>(dosage.sparse_value_ct) *
            (sizeof(uint32_t) + sizeof(uint16_t));
        if (copied_sparse_bytes + value_bytes <=
            kMaximumCopiedSparseBytes) {
          variant_states[local_variant_idx] = {
              TileVariantStorage::kSparse,
              static_cast<uint32_t>(sparse_dosages.size())};
          StoredSparseDosage stored;
          stored.common = dosage.common;
          stored.mean = dosage.mean;
          if (dosage.sparse_value_ct) {
            stored.sample_ids.assign(dosage.sparse_sample_ids,
                                     dosage.sparse_sample_ids +
                                         dosage.sparse_value_ct);
            stored.dosage16.assign(dosage.sparse_dosage16,
                                   dosage.sparse_dosage16 +
                                       dosage.sparse_value_ct);
          }
          sparse_dosages.push_back(std::move(stored));
          copied_sparse_bytes += value_bytes;
        } else {
          variant_states[local_variant_idx] = {
              TileVariantStorage::kDense,
              static_cast<uint32_t>(dense_dosages.size() /
                                    matrix->column_ct())};
          const size_t begin = dense_dosages.size();
          dense_dosages.resize(begin + matrix->column_ct(), dosage.common);
          for (uint32_t value_idx = 0; value_idx < dosage.sparse_value_ct;
               ++value_idx) {
            dense_dosages[begin + dosage.sparse_sample_ids[value_idx]] =
                dosage.sparse_dosage16[value_idx] == UINT16_MAX
                    ? dosage.mean
                    : static_cast<double>(dosage.sparse_dosage16[value_idx]) /
                          16384.0;
          }
          ++stats.densified_sparse_variant_ct;
        }
      } else {
        ++stats.dense_variant_ct;
        variant_states[local_variant_idx] = {
            TileVariantStorage::kDense,
            static_cast<uint32_t>(dense_dosages.size() /
                                  matrix->column_ct())};
        const size_t begin = dense_dosages.size();
        dense_dosages.resize(begin + matrix->column_ct());
        std::copy(dosage.dense_values,
                  dosage.dense_values + matrix->column_ct(),
                  dense_dosages.begin() + begin);
      }
    }

    std::vector<ScoreRowTask> tasks;
    uint64_t tile_edge_ct = 0;
    for (uint32_t fragment_idx = 0; fragment_idx < tiles.size();
         ++fragment_idx) {
      const auto& score_map = score_maps[fragment_idx];
      for (const auto& row : tiles[fragment_idx].rows()) {
        if (row.local_score_idx() >= score_map.size()) {
          throw std::runtime_error(
              "score-major fragment row has invalid local score index");
        }
        tasks.push_back({&row, score_map[row.local_score_idx()]});
        tile_edge_ct += row.edge_ct();
      }
    }
    std::vector<bool> score_seen(catalog->scores.size(), false);
    for (const auto& task : tasks) {
      if (score_seen[task.output_score_idx]) {
        throw std::runtime_error(
            "score-major tile contains the same output score twice");
      }
      score_seen[task.output_score_idx] = true;
    }
    const auto scoring_started = std::chrono::steady_clock::now();
    const ScoreRowWorkStats work = scoring_workers.Dispatch(
        tasks, variant_states, sparse_dosages);
    const auto plan_started = std::chrono::steady_clock::now();
    const uint32_t dense_variant_ct = static_cast<uint32_t>(
        dense_dosages.size() / matrix->column_ct());
    const BlockedDensePlan dense_plan = BuildBlockedDensePlan(
        scoring_workers.dense_edges(), dense_variant_ct, catalog, matrix);
    const auto plan_finished = std::chrono::steady_clock::now();
    if (dense_plan.weights.size() != work.dense_edge_ct) {
      throw std::runtime_error("blocked-dense plan lost score weights");
    }
    const uint32_t dense_sample_block_ct =
        dense_workers.Dispatch(dense_plan, dense_dosages, dense_variant_ct);
    const auto scoring_finished = std::chrono::steady_clock::now();
    stats.score_major_scoring_nanoseconds += static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            scoring_finished - scoring_started)
            .count());
    stats.blocked_dense_plan_nanoseconds += static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(plan_finished -
                                                             plan_started)
            .count());
    stats.blocked_dense_scoring_nanoseconds += static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            scoring_finished - plan_finished)
            .count());
    if (!dense_plan.weights.empty()) {
      ++stats.blocked_dense_tile_ct;
      stats.blocked_dense_sample_block_ct += dense_sample_block_ct;
      stats.blocked_dense_row_ct += dense_plan.rows.size();
    }
    ++stats.score_major_tile_ct;
    stats.score_major_row_ct += work.row_ct;
    stats.score_major_maximum_rows_per_tile = std::max<uint64_t>(
        stats.score_major_maximum_rows_per_tile, tasks.size());
    stats.score_major_maximum_edges_per_tile = std::max(
        stats.score_major_maximum_edges_per_tile, tile_edge_ct);
    stats.edge_ct += work.scored_edge_ct;
    stats.omitted_frequency_edge_ct += work.omitted_frequency_edge_ct;
    stats.dense_edge_ct += work.dense_edge_ct;
    stats.sparse_edge_ct += work.sparse_edge_ct;
    stats.dense_update_ct += work.dense_update_ct;
    stats.sparse_update_ct += work.sparse_update_ct;
    stats.copied_sparse_genotype_bytes += copied_sparse_bytes;
    stats.maximum_genotype_buffer_bytes = std::max(
        stats.maximum_genotype_buffer_bytes,
        static_cast<uint64_t>(dense_dosages.size()) * sizeof(double) +
            copied_sparse_bytes);
    if (progress) {
      progress->MaybeEvent(
          "score", "score_major_tiles",
          {{"tiles_processed", tile_idx + 1},
           {"tiles_total", tile_ct},
           {"referenced_variants_in_tile", referenced_in_tile},
           {"variant_groups_processed", variant_groups_processed},
           {"genotype_decodes", stats.variant_ct},
           {"weight_edges_processed", stats.edge_ct},
           {"sparse_variants", stats.sparse_variant_ct},
           {"dense_variants", stats.dense_variant_ct},
           {"sparse_weight_edges", stats.sparse_edge_ct},
           {"dense_weight_edges", stats.dense_edge_ct},
           {"sparse_score_updates", stats.sparse_update_ct},
           {"dense_score_updates", stats.dense_update_ct},
           {"score_major_rows", stats.score_major_row_ct},
           {"score_major_scoring_nanoseconds",
            stats.score_major_scoring_nanoseconds},
           {"blocked_dense_tiles", stats.blocked_dense_tile_ct},
           {"blocked_dense_sample_blocks",
            stats.blocked_dense_sample_block_ct},
           {"blocked_dense_rows", stats.blocked_dense_row_ct},
           {"blocked_dense_plan_nanoseconds",
            stats.blocked_dense_plan_nanoseconds},
           {"blocked_dense_scoring_nanoseconds",
            stats.blocked_dense_scoring_nanoseconds},
           {"maximum_genotype_buffer_bytes",
            stats.maximum_genotype_buffer_bytes},
           {"imputed_values", stats.imputed_value_ct}});
    }
    if ((tile_idx + 1) % 100 == 0) {
      std::cerr << "processed " << tile_idx + 1 << "/" << tile_ct
                << " score-major tiles\n";
    }
  }
  for (uint32_t score_idx = 0; score_idx < catalog->scores.size(); ++score_idx) {
    double* row = matrix->Row(score_idx);
    const double baseline = baselines[score_idx];
    for (uint32_t sample_idx = 0; sample_idx < matrix->column_ct(); ++sample_idx) {
      row[sample_idx] += baseline;
    }
  }
  matrix->Flush();
  if (progress) {
    progress->Event(
        "score", "fragment_scoring_complete",
        {{"fragments", fragments.size()},
         {"score_major_tiles", stats.score_major_tile_ct},
         {"variant_groups_processed", variant_groups_processed},
         {"genotype_decodes", stats.variant_ct},
         {"weight_edges_processed", stats.edge_ct},
         {"sparse_variants", stats.sparse_variant_ct},
         {"dense_variants", stats.dense_variant_ct},
         {"sparse_weight_edges", stats.sparse_edge_ct},
         {"dense_weight_edges", stats.dense_edge_ct},
         {"sparse_score_updates", stats.sparse_update_ct},
         {"dense_score_updates", stats.dense_update_ct},
         {"score_major_rows", stats.score_major_row_ct},
         {"score_major_scoring_nanoseconds",
          stats.score_major_scoring_nanoseconds},
         {"blocked_dense_tiles", stats.blocked_dense_tile_ct},
         {"blocked_dense_sample_blocks",
          stats.blocked_dense_sample_block_ct},
         {"blocked_dense_rows", stats.blocked_dense_row_ct},
         {"blocked_dense_plan_nanoseconds",
          stats.blocked_dense_plan_nanoseconds},
         {"blocked_dense_scoring_nanoseconds",
          stats.blocked_dense_scoring_nanoseconds},
         {"blocked_dense_sample_block_size",
          stats.blocked_dense_sample_block_size},
         {"copied_sparse_genotype_bytes",
          stats.copied_sparse_genotype_bytes},
         {"maximum_genotype_buffer_bytes",
          stats.maximum_genotype_buffer_bytes},
         {"densified_sparse_variants", stats.densified_sparse_variant_ct},
         {"scoring_threads", thread_ct},
         {"imputed_values", stats.imputed_value_ct}});
  }
  return stats;
}

}  // namespace pgensparsescore
