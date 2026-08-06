// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

#include "progress.h"

namespace pgensparsescore {

struct VariantIndexBuildOptions {
  std::string input_path;
  std::string source_id_column = "SOURCE_ID";
  std::string target_id_column = "TARGET_ID";
  std::string ref_column = "REF";
  std::string alt_column = "ALT";
  uint32_t block_size = 100000;
  std::string output_path;
};

void BuildVariantIndex(const VariantIndexBuildOptions& options,
                       ProgressReporter* progress = nullptr);

class VariantIndex {
 public:
  explicit VariantIndex(const std::string& path);
  ~VariantIndex();

  VariantIndex(const VariantIndex&) = delete;
  VariantIndex& operator=(const VariantIndex&) = delete;

  uint64_t variant_ct() const;
  uint64_t alias_ct() const;
  uint32_t block_size() const;
  uint32_t block_ct() const;
  uint64_t signature_lo() const;
  uint64_t signature_hi() const;
  std::optional<uint32_t> Lookup(std::string_view id) const;
  std::string_view ref(uint32_t ordinal) const;
  std::string_view alt(uint32_t ordinal) const;

 private:
  int fd_ = -1;
  void* mapping_ = nullptr;
  uint64_t mapped_bytes_ = 0;
};

}  // namespace pgensparsescore
