// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "frequency.h"
#include "mapped_matrix.h"
#include "pgen_reader.h"
#include "progress.h"
#include "score_fragment.h"
#include "types.h"
#include "variant_index.h"

namespace pgensparsescore {

std::vector<std::string> ReadScoreFragmentList(const std::string& path);

struct IndexedVariantLocation {
  uint32_t input_idx = std::numeric_limits<uint32_t>::max();
  uint32_t pgen_variant_idx = 0;

  bool present() const {
    return input_idx != std::numeric_limits<uint32_t>::max();
  }
};

struct IndexedPvarStats {
  uint32_t row_ct = 0;
  uint64_t matched_variant_ct = 0;
};

IndexedPvarStats AddIndexedPvar(
    const std::string& path, uint32_t input_idx, const VariantIndex& index,
    std::vector<IndexedVariantLocation>* locations,
    ProgressReporter* progress = nullptr);

struct LoadedScoreFragments {
  std::vector<ScoreFragmentReader> fragments;
  std::vector<std::vector<uint32_t>> score_maps;
  Catalog catalog;
  uint64_t weight_ct = 0;
  uint64_t file_byte_ct = 0;
};

LoadedScoreFragments LoadScoreFragments(
    const std::vector<std::string>& paths, const VariantIndex& index,
    const std::string& score_schema_path = "",
    ProgressReporter* progress = nullptr);

ScoreRunStats ScoreFragments(
    const VariantIndex& index,
    const std::vector<ScoreFragmentReader>& fragments,
    const std::vector<std::vector<uint32_t>>& score_maps,
    const std::vector<IndexedVariantLocation>& locations,
    const IndexedFrequencyTable* frequencies,
    MissingFrequencyPolicy missing_frequency_policy,
    const std::vector<PgenDosageReader*>& readers, Catalog* catalog,
    MappedMatrix* matrix, uint32_t thread_ct = 1,
    ProgressReporter* progress = nullptr);

}  // namespace pgensparsescore
