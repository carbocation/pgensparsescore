// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <string>
#include <unordered_map>

namespace pgensparsescore {

struct AlleleFrequency {
  std::string ref;
  std::string alt;
  double alt_dosage_mean;
};

using FrequencyTable = std::unordered_map<std::string, AlleleFrequency>;

FrequencyTable ReadFrequencyTable(const std::string& path);

}  // namespace pgensparsescore
