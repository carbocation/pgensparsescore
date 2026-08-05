// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "pgen_rans_pgenlib.h"

namespace pgensparsescore {

struct DosageView {
  bool sparse = false;
  double common = 0.0;
  double mean = 0.0;
  uint32_t missing_ct = 0;
  const double* dense_values = nullptr;
  const uint32_t* sparse_sample_ids = nullptr;
  const uint16_t* sparse_dosage16 = nullptr;
  uint32_t sparse_value_ct = 0;
};

class PgenDosageReader {
 public:
  explicit PgenDosageReader(const std::string& path);
  ~PgenDosageReader();

  PgenDosageReader(const PgenDosageReader&) = delete;
  PgenDosageReader& operator=(const PgenDosageReader&) = delete;

  uint32_t sample_ct() const { return sample_ct_; }
  uint32_t variant_ct() const { return variant_ct_; }
  const char* storage_mode_name() const;
  DosageView Read(uint32_t variant_idx,
                  std::optional<double> imputation_mean = std::nullopt);

 private:
  pgen_rans::UnifiedPgenReader reader_;
  uintptr_t* genovec_ = nullptr;
  uintptr_t* dosage_present_ = nullptr;
  uint16_t* dosage_main_ = nullptr;
  uint32_t* difflist_sample_ids_ = nullptr;
  uint32_t sample_ct_ = 0;
  uint32_t variant_ct_ = 0;
  uint32_t max_sparse_dosage_ct_ = 0;
  plink2::PgrSampleSubsetIndex pssi_{};
  std::vector<double> dense_values_;
};

}  // namespace pgensparsescore
