// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <string>
#include <vector>

#include "types.h"

namespace pgensparsescore {

std::vector<Variant> ReadPvar(const std::string& path);
std::vector<Sample> ReadPsam(const std::string& path);
Catalog CompileCatalog(const std::string& manifest_path,
                       const std::vector<Variant>& variants);

}  // namespace pgensparsescore
