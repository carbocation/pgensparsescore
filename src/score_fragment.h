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
  std::string temporary_directory;
  std::string output_path;
};

struct FragmentScore {
  std::string score_id;
  ScoreInfo info;
};

struct ScoreFragmentSummary {
  uint64_t variant_index_variant_ct = 0;
  uint32_t block_size = 0;
  uint32_t block_ct = 0;
  uint32_t score_ct = 0;
  uint64_t weight_ct = 0;
  uint64_t input_weight_ct = 0;
  uint64_t zero_weight_ct = 0;
  uint64_t excluded_weight_ct = 0;
  uint64_t duplicate_weight_ct = 0;
  uint64_t output_bytes = 0;
};

ScoreFragmentSummary CompileScoreFragment(
    const ScoreFragmentCompileOptions& options,
    ProgressReporter* progress = nullptr);

struct IndexedVariantEdges {
  uint32_t ordinal = 0;
  std::vector<Edge> edges;
};

class ScoreFragmentBlockCursor {
 public:
  ~ScoreFragmentBlockCursor();
  ScoreFragmentBlockCursor(ScoreFragmentBlockCursor&&) noexcept;
  ScoreFragmentBlockCursor& operator=(ScoreFragmentBlockCursor&&) noexcept;

  bool done() const;
  uint32_t ordinal() const;
  void AppendEdges(std::vector<Edge>* edges) const;
  void Next();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  ScoreFragmentBlockCursor();
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
  uint32_t block_size() const;
  uint32_t block_ct() const;
  uint64_t signature_lo() const;
  uint64_t signature_hi() const;
  uint64_t weight_ct() const;
  uint64_t file_bytes() const;
  const std::vector<FragmentScore>& scores() const;
  ScoreFragmentBlockCursor OpenBlock(uint32_t block_idx) const;
  void ReadBlock(uint32_t block_idx,
                 std::vector<IndexedVariantEdges>* variants) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pgensparsescore
