// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "types.h"

namespace pgensparsescore {

using VariantMap = std::unordered_map<std::string, std::string>;

std::vector<Variant> ReadPvar(const std::string& path);
std::vector<Sample> ReadPsam(const std::string& path);
VariantMap ReadVariantMap(const std::string& path);
Catalog CompileCatalog(const std::string& manifest_path,
                       const std::vector<Variant>& variants,
                       const VariantMap* variant_map = nullptr);

}  // namespace pgensparsescore
