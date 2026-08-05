// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>

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

}  // namespace pgensparsescore
