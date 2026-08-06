// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace pgensparsescore {

struct Variant {
  std::string chrom;
  std::string id;
  std::string ref;
  std::string alt;
  uint32_t pgen_variant_idx = 0;
};

struct Sample {
  std::optional<std::string> fid;
  std::string iid;
};

struct Edge {
  uint32_t score_idx;
  double beta_alt;
  bool ref_effect = false;
};

struct VariantEdges {
  uint32_t variant_idx;
  std::vector<Edge> edges;
};

struct ScoreInfo {
  std::string id;
  std::string path;
  uint64_t input_weight_ct = 0;
  uint64_t zero_weight_ct = 0;
  uint64_t excluded_weight_ct = 0;
  uint64_t duplicate_weight_ct = 0;
  uint64_t catalog_weight_ct = 0;
  uint64_t matched_weight_ct = 0;
  uint64_t missing_variant_ct = 0;
  uint64_t missing_frequency_ct = 0;
  uint64_t ref_effect_ct = 0;
  uint64_t alt_effect_ct = 0;
  double ref_effect_intercept = 0.0;
};

struct Catalog {
  std::vector<ScoreInfo> scores;
  std::vector<VariantEdges> variants;
  std::vector<double> intercepts;
};

struct ScoreRunStats {
  uint64_t variant_ct = 0;
  uint64_t edge_ct = 0;
  uint64_t sparse_variant_ct = 0;
  uint64_t dense_variant_ct = 0;
  uint64_t sparse_edge_ct = 0;
  uint64_t dense_edge_ct = 0;
  uint64_t sparse_value_ct = 0;
  uint64_t sparse_update_ct = 0;
  uint64_t dense_update_ct = 0;
  uint64_t parallel_variant_ct = 0;
  uint64_t parallel_update_ct = 0;
  uint64_t blocked_dense_block_ct = 0;
  uint64_t blocked_dense_variant_ct = 0;
  uint64_t blocked_dense_edge_ct = 0;
  uint32_t blocked_dense_maximum_variant_ct = 0;
  uint64_t blocked_dense_maximum_edge_ct = 0;
  uint64_t blocked_dense_scoring_nanoseconds = 0;
  uint64_t imputed_value_ct = 0;
  uint64_t external_frequency_variant_ct = 0;
  uint64_t cohort_frequency_variant_ct = 0;
  uint64_t missing_frequency_variant_ct = 0;
  uint64_t omitted_frequency_variant_ct = 0;
  uint64_t omitted_frequency_edge_ct = 0;
};

}  // namespace pgensparsescore
