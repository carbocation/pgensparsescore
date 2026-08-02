// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "pgenlib_read.h"

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
  DosageView Read(uint32_t variant_idx);

 private:
  plink2::PgenFileInfo pgfi_{};
  plink2::PgenReader pgr_{};
  unsigned char* pgfi_alloc_ = nullptr;
  unsigned char* pgr_alloc_ = nullptr;
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
