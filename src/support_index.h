// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <string>

#include "progress.h"
#include "variant_index.h"

namespace pgensparsescore {

enum class VariantSupport : uint8_t {
  kMissingVariant = 0,
  kMissingFrequency = 1,
  kUsable = 2,
};

struct SupportIndexBuildOptions {
  std::string variant_index_path;
  std::string pvar_path;
  std::string frequency_path;
  std::string output_path;
};

struct SupportIndexSummary {
  uint64_t variant_ct = 0;
  uint64_t missing_variant_ct = 0;
  uint64_t missing_frequency_ct = 0;
  uint64_t usable_variant_ct = 0;
  uint64_t pvar_row_ct = 0;
  uint64_t frequency_row_ct = 0;
  uint64_t output_bytes = 0;
};

SupportIndexSummary BuildSupportIndex(
    const SupportIndexBuildOptions& options,
    ProgressReporter* progress = nullptr);

class SupportIndex {
 public:
  explicit SupportIndex(const std::string& path);
  ~SupportIndex();

  SupportIndex(const SupportIndex&) = delete;
  SupportIndex& operator=(const SupportIndex&) = delete;

  uint64_t variant_ct() const;
  uint64_t signature_lo() const;
  uint64_t signature_hi() const;
  uint64_t missing_variant_ct() const;
  uint64_t missing_frequency_ct() const;
  uint64_t usable_variant_ct() const;
  uint64_t file_bytes() const;
  VariantSupport state(uint32_t ordinal) const;

 private:
  int fd_ = -1;
  void* mapping_ = nullptr;
  uint64_t mapped_bytes_ = 0;
};

}  // namespace pgensparsescore
