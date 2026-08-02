// SPDX-License-Identifier: GPL-3.0-only
#include "pgen_reader.h"

#include <algorithm>
#include <stdexcept>

#include "pgenlib_ffi_support.h"
#include "plink2_base.h"
#include "plink2_bits.h"

namespace pgensparsescore {

namespace {

void ThrowPgenError(const std::string& operation, plink2::PglErr error) {
  throw std::runtime_error(operation + " failed with pgenlib error " +
                           std::to_string(static_cast<uint32_t>(error)));
}

}  // namespace

PgenDosageReader::PgenDosageReader(const std::string& path) {
  using namespace plink2;
  PreinitPgfi(&pgfi_);
  PreinitPgr(&pgr_);
  PgenHeaderCtrl header_ctrl;
  uintptr_t alloc_cacheline_ct = 0;
  char error_buffer[kPglErrstrBufBlen]{};
  PglErr error = PgfiInitPhase1(path.c_str(), nullptr, UINT32_MAX, UINT32_MAX,
                                &header_ctrl, &pgfi_, &alloc_cacheline_ct,
                                error_buffer);
  if (error != kPglRetSuccess) {
    throw std::runtime_error("cannot open PGEN " + path + ": " + error_buffer);
  }
  sample_ct_ = pgfi_.raw_sample_ct;
  variant_ct_ = pgfi_.raw_variant_ct;
  if (!sample_ct_ || !variant_ct_) {
    throw std::runtime_error("PGEN must contain at least one sample and variant");
  }
  if (cachealigned_malloc(alloc_cacheline_ct * kCacheline, &pgfi_alloc_)) {
    throw std::bad_alloc();
  }
  uint32_t max_vrec_width = 0;
  error = PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, variant_ct_,
                         &max_vrec_width, &pgfi_, pgfi_alloc_,
                         &alloc_cacheline_ct, error_buffer);
  if (error != kPglRetSuccess) {
    throw std::runtime_error("cannot initialize PGEN index " + path + ": " +
                             error_buffer);
  }
  if (cachealigned_malloc(alloc_cacheline_ct * kCacheline, &pgr_alloc_)) {
    throw std::bad_alloc();
  }
  error = PgrInit(path.c_str(), max_vrec_width, &pgfi_, &pgr_, pgr_alloc_);
  if (error != kPglRetSuccess) {
    ThrowPgenError("PgrInit", error);
  }
  PgrClearSampleSubsetIndex(&pgr_, &pssi_);

  if (cachealigned_malloc(NypCtToVecCt(sample_ct_) * kBytesPerVec, &genovec_) ||
      cachealigned_malloc(BitCtToVecCt(sample_ct_) * kBytesPerVec,
                          &dosage_present_) ||
      cachealigned_malloc(RoundUpPow2(sample_ct_ * sizeof(uint16_t), kCacheline),
                          &dosage_main_) ||
      cachealigned_malloc(
          RoundUpPow2((static_cast<uint64_t>(sample_ct_) + 1) * sizeof(uint32_t),
                      kCacheline),
          &difflist_sample_ids_)) {
    throw std::bad_alloc();
  }
  if (sample_ct_ > 1) {
    const uint32_t suggested =
        3 * (sample_ct_ / static_cast<uint32_t>(kPglMaxDifflistLenDivisor));
    max_sparse_dosage_ct_ =
        std::min(sample_ct_ - 1, std::max<uint32_t>(1, suggested));
  }
  dense_values_.resize(sample_ct_);
}

PgenDosageReader::~PgenDosageReader() {
  using namespace plink2;
  PglErr cleanup_error = kPglRetSuccess;
  CleanupPgr(&pgr_, &cleanup_error);
  CleanupPgfi(&pgfi_, &cleanup_error);
  if (pgfi_alloc_) aligned_free(pgfi_alloc_);
  if (pgr_alloc_) aligned_free(pgr_alloc_);
  if (genovec_) aligned_free(genovec_);
  if (dosage_present_) aligned_free(dosage_present_);
  if (dosage_main_) aligned_free(dosage_main_);
  if (difflist_sample_ids_) aligned_free(difflist_sample_ids_);
}

DosageView PgenDosageReader::Read(uint32_t variant_idx) {
  using namespace plink2;
  if (variant_idx >= variant_ct_) {
    throw std::out_of_range("PGEN variant index out of range");
  }
  uint32_t dosage_ct = 0;
  uint16_t common_dosage16 = 1;
  PglErr error;
  if (max_sparse_dosage_ct_) {
    error = PgrGetDMaybeSparse(nullptr, pssi_, sample_ct_, variant_idx,
                               max_sparse_dosage_ct_, &pgr_, genovec_,
                               dosage_present_, dosage_main_, &dosage_ct,
                               &common_dosage16, difflist_sample_ids_);
  } else {
    error = PgrGetD(nullptr, pssi_, sample_ct_, variant_idx, &pgr_, genovec_,
                    dosage_present_, dosage_main_, &dosage_ct);
  }
  if (error != kPglRetSuccess) {
    ThrowPgenError("reading PGEN variant " + std::to_string(variant_idx), error);
  }

  DosageView result;
  if (max_sparse_dosage_ct_ && common_dosage16 != 1) {
    result.sparse = true;
    result.sparse_sample_ids = difflist_sample_ids_;
    result.sparse_dosage16 = dosage_main_;
    result.sparse_value_ct = dosage_ct;
    const bool common_missing = common_dosage16 == UINT16_MAX;
    uint64_t nonmissing_ct = common_missing ? 0 : sample_ct_ - dosage_ct;
    double sum = common_missing
                     ? 0.0
                     : static_cast<double>(sample_ct_ - dosage_ct) *
                           common_dosage16 / 16384.0;
    result.missing_ct = common_missing ? sample_ct_ - dosage_ct : 0;
    for (uint32_t idx = 0; idx < dosage_ct; ++idx) {
      if (dosage_main_[idx] == UINT16_MAX) {
        ++result.missing_ct;
      } else {
        ++nonmissing_ct;
        sum += static_cast<double>(dosage_main_[idx]) / 16384.0;
      }
    }
    if (!nonmissing_ct) {
      throw std::runtime_error("scored PGEN variant " +
                               std::to_string(variant_idx) +
                               " is missing in every sample");
    }
    result.mean = sum / static_cast<double>(nonmissing_ct);
    result.common = common_missing
                        ? result.mean
                        : static_cast<double>(common_dosage16) / 16384.0;
    return result;
  }

  if (Dosage16ToDoublesMeanimpute(genovec_, dosage_present_, dosage_main_,
                                  sample_ct_, dosage_ct,
                                  dense_values_.data())) {
    throw std::runtime_error("scored PGEN variant " +
                             std::to_string(variant_idx) +
                             " is missing in every sample");
  }
  result.dense_values = dense_values_.data();
  result.missing_ct = 0;
  for (uint32_t idx = 0; idx < sample_ct_; ++idx) {
    if (GetNyparrEntry(genovec_, idx) == 3 &&
        !IsSet(dosage_present_, idx)) {
      ++result.missing_ct;
    }
  }
  return result;
}

}  // namespace pgensparsescore
