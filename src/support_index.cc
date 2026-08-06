// SPDX-License-Identifier: GPL-3.0-only
#include "support_index.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "fragment_scorer.h"
#include "frequency.h"

namespace pgensparsescore {
namespace {

constexpr char kMagic[8] = {'P', 'G', 'S', 'S', 'S', 'U', 'P', 'P'};
constexpr uint32_t kVersion = 1;
constexpr uint32_t kHeaderBytes = 72;

void PutU32(unsigned char* output, uint32_t value) {
  for (uint32_t idx = 0; idx < 4; ++idx) {
    output[idx] = static_cast<unsigned char>(value >> (8 * idx));
  }
}

void PutU64(unsigned char* output, uint64_t value) {
  for (uint32_t idx = 0; idx < 8; ++idx) {
    output[idx] = static_cast<unsigned char>(value >> (8 * idx));
  }
}

uint32_t GetU32(const unsigned char* input) {
  uint32_t result = 0;
  for (uint32_t idx = 0; idx < 4; ++idx) {
    result |= static_cast<uint32_t>(input[idx]) << (8 * idx);
  }
  return result;
}

uint64_t GetU64(const unsigned char* input) {
  uint64_t result = 0;
  for (uint32_t idx = 0; idx < 8; ++idx) {
    result |= static_cast<uint64_t>(input[idx]) << (8 * idx);
  }
  return result;
}

uint64_t FileSize(int fd, const std::string& path) {
  struct stat status {};
  if (fstat(fd, &status) || status.st_size < 0) {
    throw std::runtime_error("cannot stat " + path + ": " +
                             std::strerror(errno));
  }
  return static_cast<uint64_t>(status.st_size);
}

class RemoveFileOnExit {
 public:
  explicit RemoveFileOnExit(std::string path) : path_(std::move(path)) {}
  ~RemoveFileOnExit() {
    if (!active_) return;
    std::error_code ignored;
    std::filesystem::remove(path_, ignored);
  }
  void Release() { active_ = false; }

 private:
  std::string path_;
  bool active_ = true;
};

}  // namespace

SupportIndexSummary BuildSupportIndex(const SupportIndexBuildOptions& options,
                                      ProgressReporter* progress) {
  if (std::filesystem::exists(options.output_path)) {
    throw std::runtime_error("support-index output already exists: " +
                             options.output_path);
  }
  VariantIndex index(options.variant_index_path);
  if (progress) {
    progress->Event("support_index", "start",
                    {{"variant_index_variants", index.variant_ct()}},
                    {{"variant_index", options.variant_index_path},
                     {"pvar", options.pvar_path},
                     {"frequencies", options.frequency_path},
                     {"output", options.output_path}});
  }
  std::vector<IndexedVariantLocation> locations(index.variant_ct());
  const auto pvar =
      AddIndexedPvar(options.pvar_path, 0, index, &locations, progress);
  const auto frequencies = ReadIndexedFrequencyTable(
      options.frequency_path, index, progress);
  std::vector<uint8_t> states(index.variant_ct());
  SupportIndexSummary summary;
  summary.variant_ct = index.variant_ct();
  summary.pvar_row_ct = pvar.row_ct;
  summary.frequency_row_ct = frequencies.matched_row_ct;
  for (uint64_t ordinal = 0; ordinal < index.variant_ct(); ++ordinal) {
    VariantSupport state = VariantSupport::kUsable;
    if (!locations[ordinal].present()) {
      state = VariantSupport::kMissingVariant;
      ++summary.missing_variant_ct;
    } else if (std::isnan(frequencies.alt_dosage_means[ordinal])) {
      state = VariantSupport::kMissingFrequency;
      ++summary.missing_frequency_ct;
    } else {
      ++summary.usable_variant_ct;
    }
    states[ordinal] = static_cast<uint8_t>(state);
  }
  summary.output_bytes = kHeaderBytes + summary.variant_ct;
  std::vector<unsigned char> header(kHeaderBytes, 0);
  std::memcpy(header.data(), kMagic, sizeof(kMagic));
  PutU32(header.data() + 8, kVersion);
  PutU32(header.data() + 12, kHeaderBytes);
  PutU64(header.data() + 16, summary.variant_ct);
  PutU64(header.data() + 24, index.signature_lo());
  PutU64(header.data() + 32, index.signature_hi());
  PutU64(header.data() + 40, summary.missing_variant_ct);
  PutU64(header.data() + 48, summary.missing_frequency_ct);
  PutU64(header.data() + 56, summary.usable_variant_ct);
  PutU64(header.data() + 64, summary.output_bytes);

  const std::string temporary_path = options.output_path + ".tmp";
  if (std::filesystem::exists(temporary_path)) {
    throw std::runtime_error("temporary support-index output already exists: " +
                             temporary_path);
  }
  RemoveFileOnExit cleanup(temporary_path);
  std::ofstream output(temporary_path, std::ios::binary | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot create " + temporary_path);
  output.write(reinterpret_cast<const char*>(header.data()), header.size());
  output.write(reinterpret_cast<const char*>(states.data()), states.size());
  output.close();
  if (!output) throw std::runtime_error("cannot finish " + temporary_path);
  std::filesystem::rename(temporary_path, options.output_path);
  cleanup.Release();
  if (progress) {
    progress->Event(
        "support_index", "complete",
        {{"variants", summary.variant_ct},
         {"missing_variants", summary.missing_variant_ct},
         {"missing_frequencies", summary.missing_frequency_ct},
         {"usable_variants", summary.usable_variant_ct},
         {"pvar_rows", summary.pvar_row_ct},
         {"frequency_rows", summary.frequency_row_ct},
         {"output_bytes", summary.output_bytes}});
  }
  return summary;
}

SupportIndex::SupportIndex(const std::string& path) {
  fd_ = open(path.c_str(), O_RDONLY);
  if (fd_ < 0) {
    throw std::runtime_error("cannot open support index " + path + ": " +
                             std::strerror(errno));
  }
  mapped_bytes_ = FileSize(fd_, path);
  if (mapped_bytes_ < kHeaderBytes) {
    throw std::runtime_error(path + " is not a complete support index");
  }
  mapping_ = mmap(nullptr, static_cast<size_t>(mapped_bytes_), PROT_READ,
                  MAP_SHARED, fd_, 0);
  if (mapping_ == MAP_FAILED) {
    mapping_ = nullptr;
    throw std::runtime_error("cannot map support index " + path + ": " +
                             std::strerror(errno));
  }
  const auto* bytes = static_cast<const unsigned char*>(mapping_);
  const uint64_t variant_ct = GetU64(bytes + 16);
  const uint64_t missing_variant_ct = GetU64(bytes + 40);
  const uint64_t missing_frequency_ct = GetU64(bytes + 48);
  const uint64_t usable_variant_ct = GetU64(bytes + 56);
  const uint64_t stated_bytes = GetU64(bytes + 64);
  if (std::memcmp(bytes, kMagic, sizeof(kMagic)) != 0 ||
      GetU32(bytes + 8) != kVersion || GetU32(bytes + 12) != kHeaderBytes ||
      !variant_ct || variant_ct > UINT32_MAX ||
      missing_variant_ct + missing_frequency_ct + usable_variant_ct !=
          variant_ct ||
      stated_bytes != mapped_bytes_ ||
      mapped_bytes_ != kHeaderBytes + variant_ct) {
    throw std::runtime_error(path + " has an invalid support-index header");
  }
  const auto* states = bytes + kHeaderBytes;
  for (uint64_t ordinal = 0; ordinal < variant_ct; ++ordinal) {
    if (states[ordinal] > static_cast<uint8_t>(VariantSupport::kUsable)) {
      throw std::runtime_error(path + " has an invalid support state");
    }
  }
}

SupportIndex::~SupportIndex() {
  if (mapping_) munmap(mapping_, static_cast<size_t>(mapped_bytes_));
  if (fd_ >= 0) close(fd_);
}

uint64_t SupportIndex::variant_ct() const {
  return GetU64(static_cast<const unsigned char*>(mapping_) + 16);
}
uint64_t SupportIndex::signature_lo() const {
  return GetU64(static_cast<const unsigned char*>(mapping_) + 24);
}
uint64_t SupportIndex::signature_hi() const {
  return GetU64(static_cast<const unsigned char*>(mapping_) + 32);
}
uint64_t SupportIndex::missing_variant_ct() const {
  return GetU64(static_cast<const unsigned char*>(mapping_) + 40);
}
uint64_t SupportIndex::missing_frequency_ct() const {
  return GetU64(static_cast<const unsigned char*>(mapping_) + 48);
}
uint64_t SupportIndex::usable_variant_ct() const {
  return GetU64(static_cast<const unsigned char*>(mapping_) + 56);
}
uint64_t SupportIndex::file_bytes() const { return mapped_bytes_; }

VariantSupport SupportIndex::state(uint32_t ordinal) const {
  if (ordinal >= variant_ct()) throw std::out_of_range("support-index variant");
  const auto* bytes = static_cast<const unsigned char*>(mapping_);
  return static_cast<VariantSupport>(bytes[kHeaderBytes + ordinal]);
}

}  // namespace pgensparsescore
