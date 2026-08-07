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
#include <unordered_set>
#include <vector>

#include "fragment_scorer.h"
#include "frequency.h"
#include "io.h"

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

std::vector<std::string> ReadPvarList(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) {
    throw std::runtime_error(path + " is empty");
  }
  const auto header = SplitTabs(line);
  size_t pvar_idx = header.size();
  for (size_t idx = 0; idx < header.size(); ++idx) {
    std::string name = header[idx];
    if (!name.empty() && name.front() == '#') name.erase(0, 1);
    if (name == "PVAR") {
      if (pvar_idx != header.size()) {
        throw std::runtime_error(path + " has duplicate PVAR columns");
      }
      pvar_idx = idx;
    }
  }
  if (pvar_idx == header.size()) {
    throw std::runtime_error(path + " is missing column PVAR");
  }
  const auto base = std::filesystem::absolute(path).parent_path();
  std::vector<std::string> result;
  std::unordered_set<std::string> seen;
  uint64_t line_number = 1;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty()) continue;
    const auto fields = SplitTabs(line);
    if (fields.size() <= pvar_idx || fields[pvar_idx].empty()) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has no PVAR path");
    }
    std::filesystem::path resolved(fields[pvar_idx]);
    if (resolved.is_relative()) resolved = base / resolved;
    const std::string normalized = resolved.lexically_normal().string();
    if (!seen.insert(normalized).second) {
      throw std::runtime_error(path + " contains duplicate PVAR " +
                               fields[pvar_idx]);
    }
    result.push_back(normalized);
  }
  if (result.empty()) throw std::runtime_error(path + " has no PVAR inputs");
  return result;
}

}  // namespace

SupportIndexSummary BuildSupportIndex(const SupportIndexBuildOptions& options,
                                      ProgressReporter* progress) {
  if (std::filesystem::exists(options.output_path)) {
    throw std::runtime_error("support-index output already exists: " +
                             options.output_path);
  }
  VariantIndex index(options.variant_index_path);
  std::vector<std::string> pvar_paths;
  if (!options.pvar_path.empty()) {
    pvar_paths.push_back(options.pvar_path);
  } else if (!options.pvar_list_path.empty()) {
    pvar_paths = ReadPvarList(options.pvar_list_path);
  } else {
    throw std::runtime_error("support index has no PVAR input");
  }
  if (progress) {
    progress->Event("support_index", "start",
                    {{"variant_index_variants", index.variant_ct()},
                     {"pvar_inputs", pvar_paths.size()}},
                    {{"variant_index", options.variant_index_path},
                     {"pvar", options.pvar_path},
                     {"pvar_list", options.pvar_list_path},
                     {"frequencies", options.frequency_path},
                     {"output", options.output_path}});
  }
  std::vector<IndexedVariantLocation> locations(index.variant_ct());
  SupportIndexSummary summary;
  summary.variant_ct = index.variant_ct();
  summary.pvar_input_ct = pvar_paths.size();
  for (uint32_t input_idx = 0; input_idx < pvar_paths.size(); ++input_idx) {
    const auto pvar = AddIndexedPvar(pvar_paths[input_idx], input_idx, index,
                                     &locations, progress);
    summary.pvar_row_ct += pvar.row_ct;
  }
  const auto frequencies = ReadIndexedFrequencyTable(
      options.frequency_path, index, progress);
  std::vector<uint8_t> states(index.variant_ct());
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
         {"pvar_inputs", summary.pvar_input_ct},
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
