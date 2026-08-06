// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "progress.h"
#include "types.h"
#include "variant_index.h"

namespace pgensparsescore {

struct ScoreFragmentCompileOptions {
  std::string manifest_path;
  std::string variant_index_path;
  std::string support_index_path;
  std::string temporary_directory;
  std::string output_path;
};

struct FragmentScore {
  std::string score_id;
  ScoreInfo info;
};

struct ScoreFragmentSummary {
  uint64_t variant_index_variant_ct = 0;
  uint32_t tile_size = 0;
  uint32_t tile_ct = 0;
  uint32_t score_ct = 0;
  uint64_t catalog_weight_ct = 0;
  uint64_t weight_ct = 0;
  uint64_t input_weight_ct = 0;
  uint64_t zero_weight_ct = 0;
  uint64_t excluded_weight_ct = 0;
  uint64_t duplicate_weight_ct = 0;
  uint64_t missing_variant_weight_ct = 0;
  uint64_t missing_frequency_weight_ct = 0;
  uint64_t output_bytes = 0;
};

ScoreFragmentSummary CompileScoreFragment(
    const ScoreFragmentCompileOptions& options,
    ProgressReporter* progress = nullptr);

struct ScoreMajorFragmentEdge {
  uint32_t local_variant_idx = 0;
  double beta_alt = 0.0;
  bool ref_effect = false;
};

class ScoreFragmentScoreRow {
 public:
  uint32_t local_score_idx() const { return local_score_idx_; }
  uint32_t edge_ct() const { return edge_ct_; }
  ScoreMajorFragmentEdge edge(uint32_t edge_idx) const;

 private:
  uint32_t local_score_idx_ = 0;
  uint32_t edge_ct_ = 0;
  uint32_t tile_variant_ct_ = 0;
  const unsigned char* edge_data_ = nullptr;
  const std::string* path_ = nullptr;
  friend class ScoreFragmentReader;
};

class ScoreFragmentTile {
 public:
  uint32_t tile_idx() const { return tile_idx_; }
  uint32_t first_ordinal() const { return first_ordinal_; }
  uint32_t variant_ct() const { return variant_ct_; }
  uint32_t referenced_variant_ct() const { return referenced_variant_ct_; }
  const std::vector<ScoreFragmentScoreRow>& rows() const { return rows_; }
  void OrReferencedVariants(std::vector<uint64_t>* words) const;

 private:
  uint32_t tile_idx_ = 0;
  uint32_t first_ordinal_ = 0;
  uint32_t variant_ct_ = 0;
  uint32_t referenced_variant_ct_ = 0;
  uint32_t bitmap_word_ct_ = 0;
  const unsigned char* bitmap_data_ = nullptr;
  std::vector<ScoreFragmentScoreRow> rows_;
  friend class ScoreFragmentReader;
};

class ScoreFragmentReader {
 public:
  explicit ScoreFragmentReader(const std::string& path);
  ~ScoreFragmentReader();

  ScoreFragmentReader(const ScoreFragmentReader&) = delete;
  ScoreFragmentReader& operator=(const ScoreFragmentReader&) = delete;
  ScoreFragmentReader(ScoreFragmentReader&&) noexcept;
  ScoreFragmentReader& operator=(ScoreFragmentReader&&) noexcept;

  uint64_t variant_ct() const;
  uint32_t tile_size() const;
  uint32_t tile_ct() const;
  uint64_t signature_lo() const;
  uint64_t signature_hi() const;
  uint64_t weight_ct() const;
  uint64_t file_bytes() const;
  const std::vector<FragmentScore>& scores() const;
  ScoreFragmentTile OpenTile(uint32_t tile_idx) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pgensparsescore
