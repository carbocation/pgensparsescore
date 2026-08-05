// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include "catalog.h"
#include "types.h"

namespace pgensparsescore {

struct CompiledWeight {
  uint32_t score_idx = 0;
  double weight = 0.0;
  uint8_t effect_allele_idx = 0;
};

struct CompiledVariant {
  std::string source_id;
  std::string allele0;
  std::string allele1;
  std::vector<CompiledWeight> weights;
};

struct CompiledCatalog {
  std::vector<ScoreInfo> scores;
  std::vector<CompiledVariant> variants;
  uint64_t weight_ct = 0;
};

CompiledCatalog CompileSourceCatalog(
    const std::string& manifest_path,
    const std::unordered_set<std::string>* included_source_ids = nullptr);
void WriteCompiledCatalog(const std::string& path,
                          const CompiledCatalog& catalog);
CompiledCatalog ReadCompiledCatalog(const std::string& path);
Catalog MaterializeCompiledCatalog(const CompiledCatalog& compiled,
                                   const std::vector<Variant>& variants,
                                   const VariantMap* variant_map = nullptr);

}  // namespace pgensparsescore
