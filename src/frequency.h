// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "variant_index.h"

namespace pgensparsescore {

struct AlleleFrequency {
  std::string ref;
  std::string alt;
  double alt_dosage_mean;
};

using FrequencyTable = std::unordered_map<std::string, AlleleFrequency>;

enum class MissingFrequencyPolicy { kCohort, kError, kOmit };

FrequencyTable ReadFrequencyTable(
    const std::string& path,
    const std::unordered_set<std::string>* included_ids = nullptr);

struct IndexedFrequencyTable {
  std::vector<double> alt_dosage_means;
  uint64_t input_row_ct = 0;
  uint64_t matched_row_ct = 0;
};

IndexedFrequencyTable ReadIndexedFrequencyTable(
    const std::string& path, const VariantIndex& index,
    ProgressReporter* progress = nullptr);

}  // namespace pgensparsescore
