// SPDX-License-Identifier: GPL-3.0-only
#include "variant_index.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "io.h"

namespace pgensparsescore {
namespace {

#pragma pack(push, 1)
struct IndexHeader {
  char magic[8];
  uint32_t version;
  uint32_t header_bytes;
  uint64_t variant_ct;
  uint64_t alias_ct;
  uint64_t slot_ct;
  uint32_t block_size;
  uint32_t reserved;
  uint64_t signature_lo;
  uint64_t signature_hi;
  uint64_t variant_records_offset;
  uint64_t allele_bytes_offset;
  uint64_t alias_slots_offset;
  uint64_t file_bytes;
};

struct VariantRecord {
  uint64_t allele_offset;
  uint32_t ref_bytes;
  uint32_t alt_bytes;
};

struct AliasSlot {
  uint64_t hash_lo;
  uint64_t hash_hi;
  uint32_t ordinal_plus_one;
  uint32_t reserved;
};
#pragma pack(pop)

static_assert(sizeof(IndexHeader) == 96);
static_assert(sizeof(VariantRecord) == 16);
static_assert(sizeof(AliasSlot) == 24);

constexpr char kMagic[8] = {'P', 'G', 'S', 'S', 'V', 'I', 'D', 'X'};
constexpr uint32_t kVersion = 1;
constexpr uint64_t kHashSeedLo = 1469598103934665603ULL;
constexpr uint64_t kHashSeedHi = 1099511628211ULL ^ 0x9e3779b97f4a7c15ULL;

uint64_t CheckedAdd(uint64_t lhs, uint64_t rhs, const char* description) {
  if (rhs > std::numeric_limits<uint64_t>::max() - lhs) {
    throw std::runtime_error(std::string(description) + " exceeds file limits");
  }
  return lhs + rhs;
}

uint64_t CheckedMultiply(uint64_t lhs, uint64_t rhs,
                         const char* description) {
  if (lhs && rhs > std::numeric_limits<uint64_t>::max() / lhs) {
    throw std::runtime_error(std::string(description) + " exceeds file limits");
  }
  return lhs * rhs;
}

uint64_t HashBytes(std::string_view value, uint64_t seed) {
  uint64_t hash = seed;
  for (const unsigned char byte : value) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  hash ^= hash >> 33;
  hash *= 0xff51afd7ed558ccdULL;
  hash ^= hash >> 33;
  hash *= 0xc4ceb9fe1a85ec53ULL;
  hash ^= hash >> 33;
  return hash;
}

std::pair<uint64_t, uint64_t> HashId(std::string_view value) {
  return {HashBytes(value, kHashSeedLo), HashBytes(value, kHashSeedHi)};
}

void UpdateSignature(std::string_view value, uint64_t* lo, uint64_t* hi) {
  const uint64_t length = value.size();
  for (uint32_t idx = 0; idx < 8; ++idx) {
    const unsigned char byte = static_cast<unsigned char>(length >> (8 * idx));
    *lo = (*lo ^ byte) * 1099511628211ULL;
    *hi = (*hi ^ (byte + 0x51U)) * 0x100000001b3ULL;
  }
  for (const unsigned char byte : value) {
    *lo = (*lo ^ byte) * 1099511628211ULL;
    *hi = (*hi ^ (byte + 0x9dU)) * 0x100000001b3ULL;
  }
}

using Header = std::unordered_map<std::string, size_t>;

Header ParseHeader(const std::string& line, const std::string& path) {
  Header result;
  const auto fields = SplitTabs(line);
  for (size_t idx = 0; idx < fields.size(); ++idx) {
    std::string name = fields[idx];
    if (!name.empty() && name.front() == '#') {
      name.erase(0, 1);
    }
    if (!result.emplace(name, idx).second) {
      throw std::runtime_error(path + " has duplicate column " + name);
    }
  }
  return result;
}

size_t RequireColumn(const Header& header, const std::string& name,
                     const std::string& path) {
  const auto iter = header.find(name);
  if (iter == header.end()) {
    throw std::runtime_error(path + " is missing column " + name);
  }
  return iter->second;
}

struct Columns {
  size_t source;
  size_t target;
  size_t ref;
  size_t alt;
  size_t maximum;
};

Columns ReadColumns(LineReader* reader, const VariantIndexBuildOptions& options) {
  std::string line;
  if (!reader->GetLine(&line)) {
    throw std::runtime_error(options.input_path + " is empty");
  }
  const Header header = ParseHeader(line, options.input_path);
  Columns result{RequireColumn(header, options.source_id_column,
                               options.input_path),
                 RequireColumn(header, options.target_id_column,
                               options.input_path),
                 RequireColumn(header, options.ref_column, options.input_path),
                 RequireColumn(header, options.alt_column, options.input_path),
                 0};
  result.maximum =
      std::max({result.source, result.target, result.ref, result.alt});
  return result;
}

void ValidateFields(const std::vector<std::string>& fields,
                    const Columns& columns, const std::string& path,
                    uint64_t line_number) {
  if (fields.size() <= columns.maximum) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has too few fields");
  }
  const auto& source = fields[columns.source];
  const auto& target = fields[columns.target];
  const auto& ref = fields[columns.ref];
  const auto& alt = fields[columns.alt];
  if (source.empty() || target.empty() || ref.empty() || alt.empty() ||
      ref == alt || ref.size() > UINT32_MAX || alt.size() > UINT32_MAX) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has an invalid ID or allele pair");
  }
}

uint64_t ChooseSlotCount(uint64_t alias_ct) {
  if (!alias_ct) {
    throw std::runtime_error("variant list has no variants");
  }
  const uint64_t required = CheckedAdd(
      CheckedMultiply(alias_ct, 10, "variant index hash table"), 6,
      "variant index hash table") /
                            7;
  uint64_t result = 1;
  while (result < required) {
    if (result > (uint64_t{1} << 62)) {
      throw std::runtime_error("variant index hash table is too large");
    }
    result <<= 1;
  }
  return result;
}

const IndexHeader* HeaderFrom(const void* mapping) {
  return static_cast<const IndexHeader*>(mapping);
}

const VariantRecord* RecordsFrom(const void* mapping) {
  const auto* bytes = static_cast<const unsigned char*>(mapping);
  return reinterpret_cast<const VariantRecord*>(
      bytes + HeaderFrom(mapping)->variant_records_offset);
}

const AliasSlot* SlotsFrom(const void* mapping) {
  const auto* bytes = static_cast<const unsigned char*>(mapping);
  return reinterpret_cast<const AliasSlot*>(
      bytes + HeaderFrom(mapping)->alias_slots_offset);
}

void InsertAlias(AliasSlot* slots, uint64_t slot_ct, std::string_view id,
                 uint32_t ordinal) {
  const auto [hash_lo, hash_hi] = HashId(id);
  uint64_t slot_idx = hash_lo & (slot_ct - 1);
  for (uint64_t probe = 0; probe < slot_ct; ++probe) {
    AliasSlot& slot = slots[slot_idx];
    if (!slot.ordinal_plus_one) {
      slot.hash_lo = hash_lo;
      slot.hash_hi = hash_hi;
      slot.ordinal_plus_one = ordinal + 1;
      return;
    }
    if (slot.hash_lo == hash_lo && slot.hash_hi == hash_hi) {
      if (slot.ordinal_plus_one != ordinal + 1) {
        throw std::runtime_error(
            "variant list contains a duplicate ID (or a 128-bit hash collision): " +
            std::string(id));
      }
      return;
    }
    slot_idx = (slot_idx + 1) & (slot_ct - 1);
  }
  throw std::runtime_error("variant index hash table is unexpectedly full");
}

class WritableMapping {
 public:
  WritableMapping(const std::string& path, uint64_t byte_ct) : path_(path) {
    if (byte_ct > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
      throw std::runtime_error("variant index is too large for this platform");
    }
    fd_ = open(path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0644);
    if (fd_ < 0) {
      throw std::runtime_error("cannot create " + path + ": " +
                               std::strerror(errno));
    }
    if (ftruncate(fd_, static_cast<off_t>(byte_ct))) {
      const std::string message = std::strerror(errno);
      close(fd_);
      fd_ = -1;
      std::filesystem::remove(path_);
      throw std::runtime_error("cannot size " + path + ": " + message);
    }
    mapping_ = mmap(nullptr, static_cast<size_t>(byte_ct), PROT_READ | PROT_WRITE,
                    MAP_SHARED, fd_, 0);
    if (mapping_ == MAP_FAILED) {
      mapping_ = nullptr;
      const std::string message = std::strerror(errno);
      close(fd_);
      fd_ = -1;
      std::filesystem::remove(path_);
      throw std::runtime_error("cannot map " + path + ": " + message);
    }
    byte_ct_ = byte_ct;
  }

  ~WritableMapping() {
    if (mapping_) munmap(mapping_, static_cast<size_t>(byte_ct_));
    if (fd_ >= 0) close(fd_);
    if (!keep_) {
      std::error_code ignored;
      std::filesystem::remove(path_, ignored);
    }
  }

  void* data() { return mapping_; }
  void Finish() {
    if (msync(mapping_, static_cast<size_t>(byte_ct_), MS_SYNC)) {
      throw std::runtime_error("cannot flush " + path_ + ": " +
                               std::strerror(errno));
    }
    keep_ = true;
  }

 private:
  std::string path_;
  int fd_ = -1;
  void* mapping_ = nullptr;
  uint64_t byte_ct_ = 0;
  bool keep_ = false;
};

}  // namespace

void BuildVariantIndex(const VariantIndexBuildOptions& options,
                       ProgressReporter* progress) {
  if (!options.block_size) {
    throw std::runtime_error("variant index block size must be positive");
  }
  if (std::filesystem::exists(options.output_path)) {
    throw std::runtime_error("variant index output already exists: " +
                             options.output_path);
  }

  uint64_t variant_ct = 0;
  uint64_t alias_ct = 0;
  uint64_t allele_byte_ct = 0;
  {
    LineReader reader(options.input_path);
    const Columns columns = ReadColumns(&reader, options);
    std::string line;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) continue;
      const auto fields = SplitTabs(line);
      ValidateFields(fields, columns, options.input_path, line_number);
      if (variant_ct == std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("variant index exceeds 32-bit ordinals");
      }
      ++variant_ct;
      alias_ct += fields[columns.source] == fields[columns.target] ? 1 : 2;
      allele_byte_ct = CheckedAdd(
          allele_byte_ct,
          CheckedAdd(fields[columns.ref].size(), fields[columns.alt].size(),
                     "variant index allele storage"),
          "variant index allele storage");
      if (progress && !(variant_ct % 1000000)) {
        progress->MaybeEvent("variant_index", "count",
                             {{"variants_counted", variant_ct},
                              {"aliases_counted", alias_ct},
                              {"allele_bytes", allele_byte_ct}});
      }
    }
  }
  if (!variant_ct) {
    throw std::runtime_error(options.input_path + " has no variant rows");
  }

  const uint64_t slot_ct = ChooseSlotCount(alias_ct);
  const uint64_t variant_records_offset = sizeof(IndexHeader);
  const uint64_t allele_bytes_offset = CheckedAdd(
      variant_records_offset,
      CheckedMultiply(variant_ct, sizeof(VariantRecord),
                      "variant index records"),
      "variant index records");
  const uint64_t alias_slots_offset = CheckedAdd(
      allele_bytes_offset, allele_byte_ct, "variant index allele storage");
  const uint64_t file_bytes = CheckedAdd(
      alias_slots_offset,
      CheckedMultiply(slot_ct, sizeof(AliasSlot), "variant index hash table"),
      "variant index");
  const std::string temporary = options.output_path + ".tmp";
  if (std::filesystem::exists(temporary)) {
    throw std::runtime_error("temporary variant index already exists: " +
                             temporary);
  }
  WritableMapping mapping(temporary, file_bytes);
  auto* bytes = static_cast<unsigned char*>(mapping.data());
  auto* header = reinterpret_cast<IndexHeader*>(bytes);
  std::memcpy(header->magic, kMagic, sizeof(kMagic));
  header->version = kVersion;
  header->header_bytes = sizeof(IndexHeader);
  header->variant_ct = variant_ct;
  header->alias_ct = alias_ct;
  header->slot_ct = slot_ct;
  header->block_size = options.block_size;
  header->variant_records_offset = variant_records_offset;
  header->allele_bytes_offset = allele_bytes_offset;
  header->alias_slots_offset = alias_slots_offset;
  header->file_bytes = file_bytes;
  auto* records = reinterpret_cast<VariantRecord*>(bytes + variant_records_offset);
  auto* allele_bytes = bytes + allele_bytes_offset;
  auto* slots = reinterpret_cast<AliasSlot*>(bytes + alias_slots_offset);

  uint64_t signature_lo = kHashSeedLo;
  uint64_t signature_hi = kHashSeedHi;
  uint64_t allele_offset = 0;
  uint64_t built_ct = 0;
  {
    LineReader reader(options.input_path);
    const Columns columns = ReadColumns(&reader, options);
    std::string line;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) continue;
      const auto fields = SplitTabs(line);
      ValidateFields(fields, columns, options.input_path, line_number);
      if (built_ct >= variant_ct) {
        throw std::runtime_error("variant list changed between index passes");
      }
      const auto& source = fields[columns.source];
      const auto& target = fields[columns.target];
      const auto& ref = fields[columns.ref];
      const auto& alt = fields[columns.alt];
      VariantRecord& record = records[built_ct];
      record.allele_offset = allele_offset;
      record.ref_bytes = static_cast<uint32_t>(ref.size());
      record.alt_bytes = static_cast<uint32_t>(alt.size());
      std::memcpy(allele_bytes + allele_offset, ref.data(), ref.size());
      allele_offset += ref.size();
      std::memcpy(allele_bytes + allele_offset, alt.data(), alt.size());
      allele_offset += alt.size();
      InsertAlias(slots, slot_ct, source, static_cast<uint32_t>(built_ct));
      if (target != source) {
        InsertAlias(slots, slot_ct, target, static_cast<uint32_t>(built_ct));
      }
      UpdateSignature(source, &signature_lo, &signature_hi);
      UpdateSignature(target, &signature_lo, &signature_hi);
      UpdateSignature(ref, &signature_lo, &signature_hi);
      UpdateSignature(alt, &signature_lo, &signature_hi);
      ++built_ct;
      if (progress && !(built_ct % 1000000)) {
        progress->MaybeEvent("variant_index", "build",
                             {{"variants_built", built_ct},
                              {"variants_total", variant_ct},
                              {"aliases_total", alias_ct},
                              {"output_bytes", file_bytes}});
      }
    }
  }
  if (built_ct != variant_ct || allele_offset != allele_byte_ct) {
    throw std::runtime_error("variant list changed between index passes");
  }
  header->signature_lo = signature_lo;
  header->signature_hi = signature_hi;
  mapping.Finish();
  std::filesystem::rename(temporary, options.output_path);
  if (progress) {
    progress->Event("variant_index", "complete",
                    {{"variants", variant_ct},
                     {"aliases", alias_ct},
                     {"hash_slots", slot_ct},
                     {"block_size", options.block_size},
                     {"output_bytes", file_bytes}});
  }
}

VariantIndex::VariantIndex(const std::string& path) {
  fd_ = open(path.c_str(), O_RDONLY);
  if (fd_ < 0) {
    throw std::runtime_error("cannot open variant index " + path + ": " +
                             std::strerror(errno));
  }
  struct stat status {};
  if (fstat(fd_, &status) || status.st_size < 0) {
    const std::string message = std::strerror(errno);
    close(fd_);
    fd_ = -1;
    throw std::runtime_error("cannot stat variant index " + path + ": " +
                             message);
  }
  mapped_bytes_ = static_cast<uint64_t>(status.st_size);
  if (mapped_bytes_ < sizeof(IndexHeader)) {
    close(fd_);
    fd_ = -1;
    throw std::runtime_error(path + " is not a complete variant index");
  }
  mapping_ = mmap(nullptr, static_cast<size_t>(mapped_bytes_), PROT_READ,
                  MAP_SHARED, fd_, 0);
  if (mapping_ == MAP_FAILED) {
    mapping_ = nullptr;
    const std::string message = std::strerror(errno);
    close(fd_);
    fd_ = -1;
    throw std::runtime_error("cannot map variant index " + path + ": " +
                             message);
  }
  const IndexHeader* header = HeaderFrom(mapping_);
  const bool valid_magic = std::memcmp(header->magic, kMagic, sizeof(kMagic)) == 0;
  const uint64_t expected_record_end =
      header->variant_records_offset +
      header->variant_ct * sizeof(VariantRecord);
  const bool valid = valid_magic && header->version == kVersion &&
                     header->header_bytes == sizeof(IndexHeader) &&
                     header->variant_ct && header->variant_ct <= UINT32_MAX &&
                     header->alias_ct && header->slot_ct &&
                     !(header->slot_ct & (header->slot_ct - 1)) &&
                     header->block_size &&
                     header->variant_records_offset == sizeof(IndexHeader) &&
                     expected_record_end == header->allele_bytes_offset &&
                     header->allele_bytes_offset <= header->alias_slots_offset &&
                     header->alias_slots_offset <= header->file_bytes &&
                     header->slot_ct <=
                         (header->file_bytes - header->alias_slots_offset) /
                             sizeof(AliasSlot) &&
                     header->alias_slots_offset +
                             header->slot_ct * sizeof(AliasSlot) ==
                         header->file_bytes &&
                     header->file_bytes == mapped_bytes_;
  if (!valid) {
    munmap(mapping_, static_cast<size_t>(mapped_bytes_));
    mapping_ = nullptr;
    close(fd_);
    fd_ = -1;
    throw std::runtime_error(path + " has an invalid variant-index header");
  }
  const VariantRecord* records = RecordsFrom(mapping_);
  const uint64_t allele_byte_ct =
      header->alias_slots_offset - header->allele_bytes_offset;
  for (uint64_t ordinal = 0; ordinal < header->variant_ct; ++ordinal) {
    const auto& record = records[ordinal];
    if (record.allele_offset > allele_byte_ct ||
        static_cast<uint64_t>(record.ref_bytes) + record.alt_bytes >
            allele_byte_ct - record.allele_offset ||
        !record.ref_bytes || !record.alt_bytes) {
      munmap(mapping_, static_cast<size_t>(mapped_bytes_));
      mapping_ = nullptr;
      close(fd_);
      fd_ = -1;
      throw std::runtime_error(path + " has an invalid allele record");
    }
  }
}

VariantIndex::~VariantIndex() {
  if (mapping_) munmap(mapping_, static_cast<size_t>(mapped_bytes_));
  if (fd_ >= 0) close(fd_);
}

uint64_t VariantIndex::variant_ct() const {
  return HeaderFrom(mapping_)->variant_ct;
}

uint64_t VariantIndex::alias_ct() const { return HeaderFrom(mapping_)->alias_ct; }

uint32_t VariantIndex::block_size() const {
  return HeaderFrom(mapping_)->block_size;
}

uint32_t VariantIndex::block_ct() const {
  const uint64_t count =
      (variant_ct() + block_size() - 1) / block_size();
  return static_cast<uint32_t>(count);
}

uint64_t VariantIndex::signature_lo() const {
  return HeaderFrom(mapping_)->signature_lo;
}

uint64_t VariantIndex::signature_hi() const {
  return HeaderFrom(mapping_)->signature_hi;
}

std::optional<uint32_t> VariantIndex::Lookup(std::string_view id) const {
  const IndexHeader* header = HeaderFrom(mapping_);
  const AliasSlot* slots = SlotsFrom(mapping_);
  const auto [hash_lo, hash_hi] = HashId(id);
  uint64_t slot_idx = hash_lo & (header->slot_ct - 1);
  for (uint64_t probe = 0; probe < header->slot_ct; ++probe) {
    const AliasSlot& slot = slots[slot_idx];
    if (!slot.ordinal_plus_one) return std::nullopt;
    if (slot.hash_lo == hash_lo && slot.hash_hi == hash_hi) {
      return slot.ordinal_plus_one - 1;
    }
    slot_idx = (slot_idx + 1) & (header->slot_ct - 1);
  }
  return std::nullopt;
}

std::string_view VariantIndex::ref(uint32_t ordinal) const {
  if (ordinal >= variant_ct()) throw std::out_of_range("variant ordinal");
  const IndexHeader* header = HeaderFrom(mapping_);
  const VariantRecord& record = RecordsFrom(mapping_)[ordinal];
  const auto* bytes = static_cast<const char*>(mapping_);
  return {bytes + header->allele_bytes_offset + record.allele_offset,
          record.ref_bytes};
}

std::string_view VariantIndex::alt(uint32_t ordinal) const {
  if (ordinal >= variant_ct()) throw std::out_of_range("variant ordinal");
  const IndexHeader* header = HeaderFrom(mapping_);
  const VariantRecord& record = RecordsFrom(mapping_)[ordinal];
  const auto* bytes = static_cast<const char*>(mapping_);
  return {bytes + header->allele_bytes_offset + record.allele_offset +
              record.ref_bytes,
          record.alt_bytes};
}

}  // namespace pgensparsescore
