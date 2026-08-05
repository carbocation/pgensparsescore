// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "types.h"

namespace pgensparsescore {

using VariantMap = std::unordered_map<std::string, std::string>;

struct PvarData {
  std::vector<Variant> variants;
  uint32_t row_ct = 0;
};

PvarData ReadPvar(
    const std::string& path,
    const std::unordered_set<std::string>* included_ids = nullptr);
std::vector<Sample> ReadPsam(const std::string& path);
VariantMap ReadVariantMap(const std::string& path);
Catalog CompileCatalog(const std::string& manifest_path,
                       const std::vector<Variant>& variants,
                       const VariantMap* variant_map = nullptr);

}  // namespace pgensparsescore
