// SPDX-License-Identifier: GPL-3.0-only
#include "score_fragment.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "io.h"
#include "support_index.h"

namespace pgensparsescore {
namespace {

constexpr char kMagic[8] = {'P', 'G', 'S', 'S', 'S', 'M', 'J', '3'};
constexpr uint32_t kVersion = 3;
constexpr uint32_t kHeaderBytes = 96;
constexpr char kVariantBitsMagic[8] = {'P', 'G', 'S', 'S', 'V', 'B', 'I', 'T'};
constexpr uint32_t kVariantBitsVersion = 1;
constexpr uint32_t kVariantBitsHeaderBytes = 64;
constexpr uint32_t kPreferredTileSize = 2000;
constexpr uint32_t kRefEffectMask = 0x80000000U;
constexpr uint32_t kIndexMask = 0x7fffffffU;
constexpr uint32_t kMaximumStringBytes = 1U << 28;

struct BucketWeight {
  uint32_t ordinal;
  uint32_t score_and_effect;
  double weight;
};
static_assert(sizeof(BucketWeight) == 16);

using Header = std::unordered_map<std::string, size_t>;

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

void WriteBytes(std::ostream* output, const void* bytes, size_t byte_ct,
                const std::string& path) {
  output->write(static_cast<const char*>(bytes),
                static_cast<std::streamsize>(byte_ct));
  if (!*output) throw std::runtime_error("cannot write " + path);
}

void WriteU32(std::ostream* output, uint32_t value, const std::string& path) {
  unsigned char bytes[4];
  PutU32(bytes, value);
  WriteBytes(output, bytes, sizeof(bytes), path);
}

void WriteU64(std::ostream* output, uint64_t value, const std::string& path) {
  unsigned char bytes[8];
  PutU64(bytes, value);
  WriteBytes(output, bytes, sizeof(bytes), path);
}

void WriteDouble(std::ostream* output, double value, const std::string& path) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  WriteU64(output, bits, path);
}

void WriteString(std::ostream* output, const std::string& value,
                 const std::string& path) {
  if (value.size() > UINT32_MAX) {
    throw std::runtime_error("score-fragment string is too long");
  }
  WriteU32(output, static_cast<uint32_t>(value.size()), path);
  WriteBytes(output, value.data(), value.size(), path);
}

uint64_t Position(std::ostream* output, const std::string& path) {
  const std::streampos position = output->tellp();
  if (position < 0) throw std::runtime_error("cannot seek " + path);
  return static_cast<uint64_t>(position);
}

class FileCleanup {
 public:
  explicit FileCleanup(std::string path) : path_(std::move(path)) {}
  ~FileCleanup() {
    if (!active_) return;
    std::error_code ignored;
    std::filesystem::remove(path_, ignored);
  }
  void Release() { active_ = false; }

 private:
  std::string path_;
  bool active_ = true;
};

double Fraction(double numerator, double denominator) {
  return denominator == 0.0 ? 0.0 : numerator / denominator;
}

void WriteScoreQc(const std::string& path,
                  const std::vector<FragmentScore>& scores,
                  double minimum_supported_fraction) {
  const std::string temporary = path + ".tmp";
  FileCleanup cleanup(temporary);
  std::ofstream output(temporary, std::ios::trunc);
  if (!output) throw std::runtime_error("cannot create " + temporary);
  output << "SCORE_ID\tCOLUMN_NAME\tPATH\tINPUT_WEIGHTS\tZERO_WEIGHTS"
            "\tEXCLUDED_WEIGHTS\tDUPLICATE_WEIGHTS\tCATALOG_WEIGHTS"
            "\tREFERENCE_SUPPORTED_WEIGHTS\tREFERENCE_MISSING_VARIANT_WEIGHTS"
            "\tREFERENCE_MISSING_FREQUENCY_WEIGHTS\tNONZERO_WEIGHT_L1"
            "\tNONZERO_WEIGHT_L2_SQUARED\tCATALOG_WEIGHT_L1"
            "\tCATALOG_WEIGHT_L2_SQUARED"
            "\tREFERENCE_SUPPORTED_WEIGHT_L1"
            "\tREFERENCE_SUPPORTED_WEIGHT_L2_SQUARED"
            "\tREFERENCE_SUPPORTED_FRACTION"
            "\tREFERENCE_SUPPORTED_L1_FRACTION"
            "\tREFERENCE_SUPPORTED_L2_SQUARED_FRACTION\tSTATUS\n";
  output << std::setprecision(17);
  for (const auto& score : scores) {
    const auto& info = score.info;
    const uint64_t supported = info.catalog_weight_ct -
                               info.missing_variant_ct -
                               info.missing_frequency_ct;
    const double supported_fraction =
        Fraction(static_cast<double>(supported),
                 static_cast<double>(info.catalog_weight_ct));
    output << score.score_id << '\t' << info.id << '\t' << info.path << '\t'
           << info.input_weight_ct << '\t' << info.zero_weight_ct << '\t'
           << info.excluded_weight_ct << '\t' << info.duplicate_weight_ct
           << '\t' << info.catalog_weight_ct << '\t' << supported << '\t'
           << info.missing_variant_ct << '\t' << info.missing_frequency_ct
           << '\t' << info.nonzero_weight_l1 << '\t' << info.nonzero_weight_l2
           << '\t' << info.catalog_weight_l1 << '\t'
           << info.catalog_weight_l2 << '\t' << info.supported_weight_l1
           << '\t' << info.supported_weight_l2 << '\t' << supported_fraction
           << '\t'
           << Fraction(info.supported_weight_l1, info.catalog_weight_l1)
           << '\t'
           << Fraction(info.supported_weight_l2, info.catalog_weight_l2)
           << '\t'
           << (supported_fraction >= minimum_supported_fraction
                   ? "eligible"
                   : "excluded_low_reference_coverage")
           << '\n';
  }
  output.close();
  if (!output) throw std::runtime_error("cannot finish " + temporary);
  std::filesystem::rename(temporary, path);
  cleanup.Release();
}

uint64_t WriteVariantBits(const std::string& path, const VariantIndex& index,
                          const std::vector<uint64_t>& words,
                          uint64_t referenced_variant_ct) {
  const std::string temporary = path + ".tmp";
  FileCleanup cleanup(temporary);
  std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot create " + temporary);
  WriteBytes(&output, kVariantBitsMagic, sizeof(kVariantBitsMagic), temporary);
  WriteU32(&output, kVariantBitsVersion, temporary);
  WriteU32(&output, kVariantBitsHeaderBytes, temporary);
  WriteU64(&output, index.variant_ct(), temporary);
  WriteU64(&output, index.signature_lo(), temporary);
  WriteU64(&output, index.signature_hi(), temporary);
  WriteU64(&output, referenced_variant_ct, temporary);
  const uint64_t file_bytes =
      kVariantBitsHeaderBytes + static_cast<uint64_t>(words.size()) * 8;
  WriteU64(&output, file_bytes, temporary);
  WriteU64(&output, 0, temporary);
  for (const uint64_t word : words) WriteU64(&output, word, temporary);
  output.close();
  if (!output) throw std::runtime_error("cannot finish " + temporary);
  std::filesystem::rename(temporary, path);
  cleanup.Release();
  return file_bytes;
}

Header MakeHeader(const std::vector<std::string>& fields,
                  const std::string& path) {
  Header result;
  for (size_t idx = 0; idx < fields.size(); ++idx) {
    std::string name = fields[idx];
    if (!name.empty() && name.front() == '#') name.erase(0, 1);
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

size_t FindScoreColumn(const Header& header, const std::string& path) {
  for (const char* name : {"COLUMN_NAME", "SCORE", "SCORE_ID"}) {
    const auto iter = header.find(name);
    if (iter != header.end()) return iter->second;
  }
  throw std::runtime_error(path +
                           " is missing COLUMN_NAME, SCORE, or SCORE_ID");
}

double ParseWeight(const std::string& value, const std::string& path,
                   uint64_t line_number) {
  errno = 0;
  char* end = nullptr;
  const double result = std::strtod(value.c_str(), &end);
  if (errno || end == value.c_str() || *end != '\0' ||
      !std::isfinite(result)) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has invalid finite weight: " + value);
  }
  return result;
}

class DirectoryCleanup {
 public:
  explicit DirectoryCleanup(std::filesystem::path path)
      : path_(std::move(path)) {}
  ~DirectoryCleanup() {
    std::error_code ignored;
    std::filesystem::remove_all(path_, ignored);
  }

 private:
  std::filesystem::path path_;
};

std::filesystem::path CreateTemporaryDirectory(
    const ScoreFragmentCompileOptions& options) {
  if (!options.temporary_directory.empty()) {
    const std::filesystem::path path(options.temporary_directory);
    if (!std::filesystem::create_directory(path)) {
      throw std::runtime_error("temporary directory already exists: " +
                               path.string());
    }
    return path;
  }
  const std::filesystem::path parent =
      std::filesystem::absolute(options.output_path).parent_path();
  const auto timestamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
  for (uint32_t attempt = 0; attempt < 100; ++attempt) {
    const auto path = parent /
                      (".pgss-fragment-" + std::to_string(getpid()) + "-" +
                       std::to_string(timestamp) + "-" +
                       std::to_string(attempt));
    if (std::filesystem::create_directory(path)) return path;
  }
  throw std::runtime_error("cannot create score-fragment temporary directory");
}

struct ManifestRow {
  std::string score_id;
  std::string column_name;
  std::string display_path;
  std::string resolved_path;
};

std::vector<ManifestRow> ReadManifest(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) throw std::runtime_error(path + " is empty");
  const Header header = MakeHeader(SplitTabs(line), path);
  size_t score_id_idx = static_cast<size_t>(-1);
  for (const char* name : {"SCORE_ID", "SCORE", "COLUMN_NAME"}) {
    const auto iter = header.find(name);
    if (iter != header.end()) {
      score_id_idx = iter->second;
      break;
    }
  }
  if (score_id_idx == static_cast<size_t>(-1)) {
    throw std::runtime_error(path +
                             " is missing SCORE_ID, SCORE, or COLUMN_NAME");
  }
  const size_t column_idx = FindScoreColumn(header, path);
  const size_t path_idx = RequireColumn(header, "PATH", path);
  const size_t maximum = std::max({score_id_idx, column_idx, path_idx});
  const std::filesystem::path base =
      std::filesystem::absolute(path).parent_path();
  std::unordered_set<std::string> seen_ids;
  std::unordered_set<std::string> seen_columns;
  std::vector<ManifestRow> result;
  uint64_t line_number = 1;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty()) continue;
    const auto fields = SplitTabs(line);
    if (fields.size() <= maximum) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has too few fields");
    }
    if (fields[score_id_idx].empty() || fields[column_idx].empty() ||
        !seen_ids.insert(fields[score_id_idx]).second ||
        !seen_columns.insert(fields[column_idx]).second) {
      throw std::runtime_error(path +
                               " has an empty or duplicate score ID/column");
    }
    std::filesystem::path resolved(fields[path_idx]);
    if (resolved.is_relative()) resolved = base / resolved;
    result.push_back({fields[score_id_idx], fields[column_idx],
                      fields[path_idx], resolved.lexically_normal().string()});
  }
  if (result.empty()) throw std::runtime_error(path + " has no scores");
  return result;
}

std::vector<std::string> ReadVariantBitsList(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) throw std::runtime_error(path + " is empty");
  const Header header = MakeHeader(SplitTabs(line), path);
  size_t bits_idx = static_cast<size_t>(-1);
  for (const char* name : {"VARIANT_BITS", "PATH"}) {
    const auto iter = header.find(name);
    if (iter != header.end()) {
      bits_idx = iter->second;
      break;
    }
  }
  if (bits_idx == static_cast<size_t>(-1)) {
    throw std::runtime_error(path + " is missing VARIANT_BITS or PATH");
  }
  const std::filesystem::path base =
      std::filesystem::absolute(path).parent_path();
  std::unordered_set<std::string> seen;
  std::vector<std::string> result;
  uint64_t line_number = 1;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty()) continue;
    const auto fields = SplitTabs(line);
    if (fields.size() <= bits_idx || fields[bits_idx].empty()) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has no variant-bitset path");
    }
    std::filesystem::path resolved(fields[bits_idx]);
    if (resolved.is_relative()) resolved = base / resolved;
    const std::string normalized = resolved.lexically_normal().string();
    if (!seen.insert(normalized).second) {
      throw std::runtime_error(path + " contains duplicate variant bitset " +
                               fields[bits_idx]);
    }
    result.push_back(normalized);
  }
  if (result.empty()) throw std::runtime_error(path + " has no variant bitsets");
  return result;
}

void FillHeader(unsigned char* header, const VariantIndex& index,
                uint32_t tile_size, uint32_t tile_ct, uint32_t score_ct,
                uint64_t weight_ct, uint64_t score_records_offset,
                uint64_t tile_directory_offset, uint64_t tile_data_offset,
                uint64_t file_bytes) {
  std::memset(header, 0, kHeaderBytes);
  std::memcpy(header, kMagic, sizeof(kMagic));
  PutU32(header + 8, kVersion);
  PutU32(header + 12, kHeaderBytes);
  PutU64(header + 16, index.variant_ct());
  PutU32(header + 24, tile_size);
  PutU32(header + 28, tile_ct);
  PutU64(header + 32, index.signature_lo());
  PutU64(header + 40, index.signature_hi());
  PutU32(header + 48, score_ct);
  PutU64(header + 56, weight_ct);
  PutU64(header + 64, score_records_offset);
  PutU64(header + 72, tile_directory_offset);
  PutU64(header + 80, tile_data_offset);
  PutU64(header + 88, file_bytes);
}

uint64_t ReadFileSize(int fd, const std::string& path) {
  struct stat status {};
  if (fstat(fd, &status) || status.st_size < 0) {
    throw std::runtime_error("cannot stat " + path + ": " +
                             std::strerror(errno));
  }
  return static_cast<uint64_t>(status.st_size);
}

const unsigned char* ReadString(const unsigned char* cursor,
                                const unsigned char* end,
                                const std::string& path,
                                std::string* value) {
  if (end - cursor < 4) throw std::runtime_error(path + " is truncated");
  const uint32_t byte_ct = GetU32(cursor);
  cursor += 4;
  if (byte_ct > kMaximumStringBytes ||
      static_cast<uint64_t>(end - cursor) < byte_ct) {
    throw std::runtime_error(path + " has an invalid string record");
  }
  value->assign(reinterpret_cast<const char*>(cursor), byte_ct);
  return cursor + byte_ct;
}

uint32_t Popcount(uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
  return static_cast<uint32_t>(__builtin_popcountll(value));
#else
  uint32_t result = 0;
  while (value) {
    value &= value - 1;
    ++result;
  }
  return result;
#endif
}

}  // namespace

ScoreFragmentSummary CompileScoreFragment(
    const ScoreFragmentCompileOptions& options, ProgressReporter* progress) {
  if (!std::isfinite(options.minimum_supported_fraction) ||
      options.minimum_supported_fraction < 0.0 ||
      options.minimum_supported_fraction > 1.0) {
    throw std::runtime_error(
        "minimum supported fraction must be from 0 through 1");
  }
  if (options.minimum_supported_fraction > 0.0 &&
      options.support_index_path.empty()) {
    throw std::runtime_error(
        "minimum supported fraction requires a support index");
  }
  if (std::filesystem::exists(options.output_path)) {
    throw std::runtime_error("score-fragment output already exists: " +
                             options.output_path);
  }
  for (const std::string suffix : {".score_qc.tsv", ".variants.bits"}) {
    if (std::filesystem::exists(options.output_path + suffix)) {
      throw std::runtime_error("score-fragment sidecar already exists: " +
                               options.output_path + suffix);
    }
  }
  VariantIndex index(options.variant_index_path);
  std::unique_ptr<SupportIndex> support;
  if (!options.support_index_path.empty()) {
    support = std::make_unique<SupportIndex>(options.support_index_path);
    if (support->variant_ct() != index.variant_ct() ||
        support->signature_lo() != index.signature_lo() ||
        support->signature_hi() != index.signature_hi()) {
      throw std::runtime_error(
          "support index was built for a different variant index");
    }
  }
  const auto manifest = ReadManifest(options.manifest_path);
  if (manifest.size() > kIndexMask) {
    throw std::runtime_error("score fragment contains too many scores");
  }
  const uint32_t tile_size =
      std::min(kPreferredTileSize, index.block_size());
  if (!tile_size || index.block_size() % tile_size) {
    throw std::runtime_error(
        "variant-index block size must be divisible by the scoring tile size");
  }
  const uint32_t tile_ct = static_cast<uint32_t>(
      (index.variant_ct() + tile_size - 1) / tile_size);
  if (progress) {
    progress->Event("fragment", "start",
                    {{"score_files_total", manifest.size()},
                     {"variant_index_variants", index.variant_ct()},
                     {"scoring_tile_size", tile_size},
                     {"scoring_tiles", tile_ct}},
                    {{"manifest", options.manifest_path},
                     {"variant_index", options.variant_index_path},
                     {"support_index", options.support_index_path},
                     {"output", options.output_path}});
  }

  const auto temp_directory = CreateTemporaryDirectory(options);
  DirectoryCleanup cleanup(temp_directory);
  std::vector<std::unique_ptr<std::ofstream>> buckets(index.block_ct());
  std::vector<uint64_t> bucket_weight_ct(index.block_ct(), 0);
  std::vector<uint64_t> referenced_variant_words(
      (index.variant_ct() + 63) / 64, 0);
  std::vector<FragmentScore> scores;
  scores.reserve(manifest.size());
  ScoreFragmentSummary summary;
  summary.variant_index_variant_ct = index.variant_ct();
  summary.tile_size = tile_size;
  summary.tile_ct = tile_ct;
  summary.input_score_ct = static_cast<uint32_t>(manifest.size());

  for (uint32_t manifest_idx = 0; manifest_idx < manifest.size();
       ++manifest_idx) {
    const auto& item = manifest[manifest_idx];
    ScoreInfo info;
    info.id = item.column_name;
    info.path = item.display_path;
    LineReader reader(item.resolved_path);
    std::string line;
    if (!reader.GetLine(&line)) {
      throw std::runtime_error(item.resolved_path + " is empty");
    }
    const Header header = MakeHeader(SplitTabs(line), item.resolved_path);
    const size_t snp_idx = RequireColumn(header, "SNP", item.resolved_path);
    const size_t effect_idx =
        RequireColumn(header, "EFFECT_ALLELE", item.resolved_path);
    const size_t other_idx =
        RequireColumn(header, "OTHER_ALLELE", item.resolved_path);
    const size_t weight_idx =
        RequireColumn(header, "EFFECT_ALLELE_WEIGHT", item.resolved_path);
    const size_t maximum =
        std::max({snp_idx, effect_idx, other_idx, weight_idx});
    std::unordered_set<uint32_t> seen_variants;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) continue;
      const auto fields = SplitTabs(line);
      if (fields.size() <= maximum) {
        throw std::runtime_error(item.resolved_path + ": line " +
                                 std::to_string(line_number) +
                                 " has too few fields");
      }
      ++info.input_weight_ct;
      ++summary.input_weight_ct;
      const double weight =
          ParseWeight(fields[weight_idx], item.resolved_path, line_number);
      if (weight == 0.0) {
        ++info.zero_weight_ct;
        ++summary.zero_weight_ct;
        continue;
      }
      const double absolute_weight = std::abs(weight);
      const double squared_weight = weight * weight;
      info.nonzero_weight_l1 += absolute_weight;
      info.nonzero_weight_l2 += squared_weight;
      const auto ordinal = index.Lookup(fields[snp_idx]);
      if (!ordinal) {
        ++info.excluded_weight_ct;
        ++summary.excluded_weight_ct;
        continue;
      }
      const std::string_view ref = index.ref(*ordinal);
      const std::string_view alt = index.alt(*ordinal);
      const std::string& effect = fields[effect_idx];
      const std::string& other = fields[other_idx];
      const bool alt_effect = effect == alt && other == ref;
      const bool ref_effect = effect == ref && other == alt;
      if (!alt_effect && !ref_effect) {
        throw std::runtime_error(
            item.resolved_path + ": line " + std::to_string(line_number) +
            " has alleles " + effect + "/" + other +
            " which disagree with indexed REF/ALT " + std::string(ref) + "/" +
            std::string(alt) + " for " + fields[snp_idx]);
      }
      if (!seen_variants.insert(*ordinal).second) {
        ++info.duplicate_weight_ct;
        ++summary.duplicate_weight_ct;
      }
      ++info.catalog_weight_ct;
      ++summary.catalog_weight_ct;
      info.catalog_weight_l1 += absolute_weight;
      info.catalog_weight_l2 += squared_weight;
      if (support) {
        const VariantSupport state = support->state(*ordinal);
        if (state == VariantSupport::kMissingVariant) {
          ++info.missing_variant_ct;
          ++summary.missing_variant_weight_ct;
          continue;
        }
        if (state == VariantSupport::kMissingFrequency) {
          ++info.matched_weight_ct;
          ++info.missing_frequency_ct;
          if (ref_effect) {
            ++info.ref_effect_ct;
          } else {
            ++info.alt_effect_ct;
          }
          ++summary.missing_frequency_weight_ct;
          continue;
        }
      }
      info.supported_weight_l1 += absolute_weight;
      info.supported_weight_l2 += squared_weight;
      const uint32_t block_idx = *ordinal / index.block_size();
      if (!buckets[block_idx]) {
        const auto path = temp_directory /
                          ("block-" + std::to_string(block_idx) + ".bin");
        buckets[block_idx] = std::make_unique<std::ofstream>(
            path, std::ios::binary | std::ios::trunc);
        if (!*buckets[block_idx]) {
          throw std::runtime_error("cannot create " + path.string());
        }
      }
      const BucketWeight record{
          *ordinal,
          manifest_idx | (ref_effect ? kRefEffectMask : 0U), weight};
      buckets[block_idx]->write(reinterpret_cast<const char*>(&record),
                                sizeof(record));
      if (!*buckets[block_idx]) {
        throw std::runtime_error("cannot write score-fragment bucket");
      }
      ++bucket_weight_ct[block_idx];
      ++summary.supported_weight_ct;
      if (progress && !(summary.input_weight_ct % 1000000)) {
        progress->MaybeEvent(
            "fragment", "parse_weights",
            {{"score_files_processed", manifest_idx},
             {"score_files_total", manifest.size()},
             {"input_weights", summary.input_weight_ct},
             {"supported_weights", summary.supported_weight_ct},
             {"catalog_weights", summary.catalog_weight_ct},
             {"missing_variant_weights",
              summary.missing_variant_weight_ct},
             {"missing_frequency_weights",
              summary.missing_frequency_weight_ct},
             {"zero_weights", summary.zero_weight_ct},
             {"excluded_weights", summary.excluded_weight_ct},
             {"duplicate_weights", summary.duplicate_weight_ct}},
            {{"current_score", item.score_id}});
      }
    }
    scores.push_back({item.score_id, std::move(info)});
    if (progress) {
      progress->MaybeEvent(
          "fragment", "parse_weights",
          {{"score_files_processed", manifest_idx + 1},
           {"score_files_total", manifest.size()},
           {"input_weights", summary.input_weight_ct},
           {"supported_weights", summary.supported_weight_ct},
           {"catalog_weights", summary.catalog_weight_ct},
           {"missing_variant_weights", summary.missing_variant_weight_ct},
           {"missing_frequency_weights",
            summary.missing_frequency_weight_ct},
           {"zero_weights", summary.zero_weight_ct},
           {"excluded_weights", summary.excluded_weight_ct},
           {"duplicate_weights", summary.duplicate_weight_ct}},
          {{"current_score", item.score_id}});
    }
  }
  for (auto& bucket : buckets) {
    if (!bucket) continue;
    bucket->close();
    if (!*bucket) throw std::runtime_error("cannot finish fragment bucket");
    bucket.reset();
  }

  std::vector<uint32_t> included_score_map(scores.size(), UINT32_MAX);
  std::vector<FragmentScore> included_scores;
  included_scores.reserve(scores.size());
  for (uint32_t score_idx = 0; score_idx < scores.size(); ++score_idx) {
    const auto& info = scores[score_idx].info;
    const uint64_t supported = info.catalog_weight_ct -
                               info.missing_variant_ct -
                               info.missing_frequency_ct;
    const double supported_fraction =
        Fraction(static_cast<double>(supported),
                 static_cast<double>(info.catalog_weight_ct));
    if (supported_fraction < options.minimum_supported_fraction) continue;
    included_score_map[score_idx] =
        static_cast<uint32_t>(included_scores.size());
    included_scores.push_back(scores[score_idx]);
    summary.weight_ct += supported;
  }
  summary.score_ct = static_cast<uint32_t>(included_scores.size());
  summary.excluded_score_ct = summary.input_score_ct - summary.score_ct;

  const std::string temporary_output = options.output_path + ".tmp";
  if (std::filesystem::exists(temporary_output)) {
    throw std::runtime_error("temporary fragment output already exists: " +
                             temporary_output);
  }
  FileCleanup output_cleanup(temporary_output);
  std::ofstream output(temporary_output, std::ios::binary | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot create " + temporary_output);
  std::vector<unsigned char> header(kHeaderBytes, 0);
  WriteBytes(&output, header.data(), header.size(), temporary_output);
  const uint64_t score_records_offset = Position(&output, temporary_output);
  for (const auto& score : included_scores) {
    WriteString(&output, score.score_id, temporary_output);
    WriteString(&output, score.info.id, temporary_output);
    WriteString(&output, score.info.path, temporary_output);
    WriteU64(&output, score.info.input_weight_ct, temporary_output);
    WriteU64(&output, score.info.zero_weight_ct, temporary_output);
    WriteU64(&output, score.info.excluded_weight_ct, temporary_output);
    WriteU64(&output, score.info.duplicate_weight_ct, temporary_output);
    WriteU64(&output, score.info.catalog_weight_ct, temporary_output);
    WriteU64(&output, score.info.matched_weight_ct, temporary_output);
    WriteU64(&output, score.info.missing_variant_ct, temporary_output);
    WriteU64(&output, score.info.missing_frequency_ct, temporary_output);
    WriteU64(&output, score.info.alt_effect_ct, temporary_output);
    WriteU64(&output, score.info.ref_effect_ct, temporary_output);
    WriteDouble(&output, score.info.nonzero_weight_l1, temporary_output);
    WriteDouble(&output, score.info.nonzero_weight_l2, temporary_output);
    WriteDouble(&output, score.info.catalog_weight_l1, temporary_output);
    WriteDouble(&output, score.info.catalog_weight_l2, temporary_output);
    WriteDouble(&output, score.info.supported_weight_l1, temporary_output);
    WriteDouble(&output, score.info.supported_weight_l2, temporary_output);
  }
  const uint64_t tile_directory_offset = Position(&output, temporary_output);
  std::vector<uint64_t> tile_offsets(static_cast<uint64_t>(tile_ct) + 1, 0);
  for (uint32_t idx = 0; idx <= tile_ct; ++idx) {
    WriteU64(&output, 0, temporary_output);
  }
  const uint64_t tile_data_offset = Position(&output, temporary_output);
  const uint32_t tiles_per_index_block = index.block_size() / tile_size;
  uint64_t weights_written = 0;
  uint32_t next_tile_idx = 0;
  for (uint32_t block_idx = 0; block_idx < index.block_ct(); ++block_idx) {
    std::vector<BucketWeight> records(bucket_weight_ct[block_idx]);
    if (!records.empty()) {
      const auto path = temp_directory /
                        ("block-" + std::to_string(block_idx) + ".bin");
      std::ifstream input(path, std::ios::binary);
      if (!input) throw std::runtime_error("cannot open " + path.string());
      input.read(reinterpret_cast<char*>(records.data()),
                 static_cast<std::streamsize>(records.size() *
                                              sizeof(records[0])));
      if (!input || input.peek() != std::ifstream::traits_type::eof()) {
        throw std::runtime_error("invalid score-fragment bucket " +
                                 path.string());
      }
      records.erase(
          std::remove_if(records.begin(), records.end(),
                         [&included_score_map](const BucketWeight& record) {
                           const uint32_t score_idx =
                               record.score_and_effect & kIndexMask;
                           return included_score_map[score_idx] == UINT32_MAX;
                         }),
          records.end());
      for (auto& record : records) {
        const uint32_t old_score_idx = record.score_and_effect & kIndexMask;
        record.score_and_effect =
            included_score_map[old_score_idx] |
            (record.score_and_effect & kRefEffectMask);
      }
      std::stable_sort(
          records.begin(), records.end(),
          [tile_size](const BucketWeight& lhs, const BucketWeight& rhs) {
            const uint32_t lhs_tile = lhs.ordinal / tile_size;
            const uint32_t rhs_tile = rhs.ordinal / tile_size;
            if (lhs_tile != rhs_tile) return lhs_tile < rhs_tile;
            const uint32_t lhs_score = lhs.score_and_effect & kIndexMask;
            const uint32_t rhs_score = rhs.score_and_effect & kIndexMask;
            if (lhs_score != rhs_score) return lhs_score < rhs_score;
            return lhs.ordinal < rhs.ordinal;
          });
    }
    const uint32_t block_first_tile = block_idx * tiles_per_index_block;
    const uint32_t block_tile_end =
        std::min(tile_ct, block_first_tile + tiles_per_index_block);
    if (next_tile_idx != block_first_tile) {
      throw std::runtime_error("score-fragment tile alignment changed");
    }
    size_t record_begin = 0;
    for (; next_tile_idx < block_tile_end; ++next_tile_idx) {
      tile_offsets[next_tile_idx] = Position(&output, temporary_output);
      const uint32_t first_ordinal = next_tile_idx * tile_size;
      const uint32_t tile_variant_ct = static_cast<uint32_t>(
          std::min<uint64_t>(tile_size, index.variant_ct() - first_ordinal));
      size_t record_end = record_begin;
      while (record_end < records.size() &&
             records[record_end].ordinal / tile_size == next_tile_idx) {
        ++record_end;
      }
      const uint32_t bitmap_word_ct = (tile_variant_ct + 63) / 64;
      std::vector<uint64_t> bitmap(bitmap_word_ct, 0);
      uint32_t score_row_ct = 0;
      uint32_t previous_score = UINT32_MAX;
      for (size_t idx = record_begin; idx < record_end; ++idx) {
        const uint32_t local_variant_idx =
            records[idx].ordinal - first_ordinal;
        if (local_variant_idx >= tile_variant_ct) {
          throw std::runtime_error("invalid score-fragment tile contents");
        }
        bitmap[local_variant_idx / 64] |= uint64_t{1}
                                                 << (local_variant_idx % 64);
        referenced_variant_words[records[idx].ordinal / 64] |=
            uint64_t{1} << (records[idx].ordinal % 64);
        const uint32_t score_idx = records[idx].score_and_effect & kIndexMask;
        if (score_idx != previous_score) {
          ++score_row_ct;
          previous_score = score_idx;
        }
      }
      uint32_t referenced_variant_ct = 0;
      for (const uint64_t word : bitmap) {
        referenced_variant_ct += Popcount(word);
      }
      WriteU32(&output, tile_variant_ct, temporary_output);
      WriteU32(&output, referenced_variant_ct, temporary_output);
      WriteU32(&output, score_row_ct, temporary_output);
      WriteU32(&output, bitmap_word_ct, temporary_output);
      for (const uint64_t word : bitmap) {
        WriteU64(&output, word, temporary_output);
      }
      size_t score_begin = record_begin;
      while (score_begin < record_end) {
        const uint32_t score_idx =
            records[score_begin].score_and_effect & kIndexMask;
        size_t score_end = score_begin + 1;
        while (score_end < record_end &&
               (records[score_end].score_and_effect & kIndexMask) ==
                   score_idx) {
          ++score_end;
        }
        if (score_end - score_begin > UINT32_MAX) {
          throw std::runtime_error("score-fragment tile row is too large");
        }
        WriteU32(&output, score_idx, temporary_output);
        WriteU32(&output, static_cast<uint32_t>(score_end - score_begin),
                 temporary_output);
        for (size_t idx = score_begin; idx < score_end; ++idx) {
          const uint32_t local_variant_idx =
              records[idx].ordinal - first_ordinal;
          const uint32_t variant_and_effect =
              local_variant_idx |
              (records[idx].score_and_effect & kRefEffectMask);
          WriteU32(&output, variant_and_effect, temporary_output);
          WriteDouble(&output, records[idx].weight, temporary_output);
        }
        weights_written += score_end - score_begin;
        score_begin = score_end;
      }
      record_begin = record_end;
      if (progress) {
        progress->MaybeEvent(
            "fragment", "serialize_tiles",
            {{"tiles_written", next_tile_idx + 1},
             {"tiles_total", tile_ct},
             {"weights_written", weights_written},
             {"weights_total", summary.weight_ct},
             {"output_bytes", Position(&output, temporary_output)}});
      }
    }
    if (record_begin != records.size()) {
      throw std::runtime_error("score-fragment bucket crossed an index block");
    }
  }
  if (next_tile_idx != tile_ct) {
    throw std::runtime_error("score-fragment tile count changed");
  }
  tile_offsets[tile_ct] = Position(&output, temporary_output);
  if (weights_written != summary.weight_ct) {
    throw std::runtime_error("fragment weight count changed during serialization");
  }
  if (summary.catalog_weight_ct !=
      summary.supported_weight_ct + summary.missing_variant_weight_ct +
          summary.missing_frequency_weight_ct) {
    throw std::runtime_error("fragment support counts do not add up");
  }
  summary.output_bytes = tile_offsets.back();
  output.seekp(static_cast<std::streamoff>(tile_directory_offset));
  for (const uint64_t offset : tile_offsets) {
    WriteU64(&output, offset, temporary_output);
  }
  FillHeader(header.data(), index, tile_size, tile_ct, included_scores.size(),
             summary.weight_ct, score_records_offset, tile_directory_offset,
             tile_data_offset, summary.output_bytes);
  output.seekp(0);
  WriteBytes(&output, header.data(), header.size(), temporary_output);
  output.close();
  if (!output) throw std::runtime_error("cannot finish " + temporary_output);
  std::filesystem::rename(temporary_output, options.output_path);
  output_cleanup.Release();
  for (const uint64_t word : referenced_variant_words) {
    summary.referenced_variant_ct += Popcount(word);
  }
  WriteScoreQc(options.output_path + ".score_qc.tsv", scores,
               options.minimum_supported_fraction);
  summary.variant_bitset_bytes = WriteVariantBits(
      options.output_path + ".variants.bits", index,
      referenced_variant_words, summary.referenced_variant_ct);
  if (progress) {
    progress->Event("fragment", "complete",
                    {{"input_scores", summary.input_score_ct},
                     {"scores", summary.score_ct},
                     {"excluded_scores", summary.excluded_score_ct},
                     {"input_weights", summary.input_weight_ct},
                     {"supported_weights", summary.supported_weight_ct},
                     {"retained_weights", summary.weight_ct},
                     {"catalog_weights", summary.catalog_weight_ct},
                     {"missing_variant_weights",
                      summary.missing_variant_weight_ct},
                     {"missing_frequency_weights",
                      summary.missing_frequency_weight_ct},
                     {"zero_weights", summary.zero_weight_ct},
                     {"excluded_weights", summary.excluded_weight_ct},
                     {"duplicate_weights", summary.duplicate_weight_ct},
                     {"scoring_tile_size", summary.tile_size},
                     {"scoring_tiles", summary.tile_ct},
                     {"output_bytes", summary.output_bytes}});
  }
  return summary;
}

VariantBitsMergeSummary MergeVariantBits(const std::string& list_path,
                                         const std::string& output_path) {
  if (std::filesystem::exists(output_path)) {
    throw std::runtime_error("variant-bitset output already exists: " +
                             output_path);
  }
  const auto paths = ReadVariantBitsList(list_path);
  VariantBitsMergeSummary summary;
  summary.input_ct = paths.size();
  uint64_t signature_lo = 0;
  uint64_t signature_hi = 0;
  std::vector<uint64_t> merged;
  for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx) {
    const std::string& path = paths[path_idx];
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open " + path);
    unsigned char header[kVariantBitsHeaderBytes];
    input.read(reinterpret_cast<char*>(header), sizeof(header));
    if (!input || std::memcmp(header, kVariantBitsMagic,
                              sizeof(kVariantBitsMagic)) ||
        GetU32(header + 8) != kVariantBitsVersion ||
        GetU32(header + 12) != kVariantBitsHeaderBytes) {
      throw std::runtime_error(path + " is not a supported variant bitset");
    }
    const uint64_t variant_ct = GetU64(header + 16);
    const uint64_t current_signature_lo = GetU64(header + 24);
    const uint64_t current_signature_hi = GetU64(header + 32);
    const uint64_t stated_referenced_ct = GetU64(header + 40);
    const uint64_t stated_file_bytes = GetU64(header + 48);
    const uint64_t word_ct = (variant_ct + 63) / 64;
    const uint64_t expected_file_bytes =
        kVariantBitsHeaderBytes + word_ct * 8;
    if (!variant_ct || stated_file_bytes != expected_file_bytes ||
        std::filesystem::file_size(path) != expected_file_bytes) {
      throw std::runtime_error(path + " has inconsistent dimensions");
    }
    if (!path_idx) {
      summary.variant_ct = variant_ct;
      signature_lo = current_signature_lo;
      signature_hi = current_signature_hi;
      merged.assign(word_ct, 0);
    } else if (variant_ct != summary.variant_ct ||
               current_signature_lo != signature_lo ||
               current_signature_hi != signature_hi) {
      throw std::runtime_error(
          path + " was built for a different variant index");
    }
    uint64_t observed_referenced_ct = 0;
    constexpr uint64_t kReadWordCt = 65536;
    std::vector<unsigned char> bytes(kReadWordCt * 8);
    for (uint64_t word_begin = 0; word_begin < word_ct;
         word_begin += kReadWordCt) {
      const uint64_t current_word_ct =
          std::min(kReadWordCt, word_ct - word_begin);
      input.read(reinterpret_cast<char*>(bytes.data()),
                 static_cast<std::streamsize>(current_word_ct * 8));
      if (!input) throw std::runtime_error(path + " is truncated");
      for (uint64_t local_idx = 0; local_idx < current_word_ct; ++local_idx) {
        const uint64_t word_idx = word_begin + local_idx;
        const uint64_t word = GetU64(bytes.data() + local_idx * 8);
        if (word_idx + 1 == word_ct && variant_ct % 64 &&
            (word >> (variant_ct % 64))) {
          throw std::runtime_error(path +
                                   " sets bits beyond its variant count");
        }
        observed_referenced_ct += Popcount(word);
        merged[word_idx] |= word;
      }
    }
    if (observed_referenced_ct != stated_referenced_ct ||
        input.peek() != std::ifstream::traits_type::eof()) {
      throw std::runtime_error(path + " has inconsistent referenced variants");
    }
  }
  for (const uint64_t word : merged) {
    summary.referenced_variant_ct += Popcount(word);
  }
  summary.output_bytes =
      kVariantBitsHeaderBytes + static_cast<uint64_t>(merged.size()) * 8;

  const std::string temporary = output_path + ".tmp";
  FileCleanup cleanup(temporary);
  std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot create " + temporary);
  WriteBytes(&output, kVariantBitsMagic, sizeof(kVariantBitsMagic), temporary);
  WriteU32(&output, kVariantBitsVersion, temporary);
  WriteU32(&output, kVariantBitsHeaderBytes, temporary);
  WriteU64(&output, summary.variant_ct, temporary);
  WriteU64(&output, signature_lo, temporary);
  WriteU64(&output, signature_hi, temporary);
  WriteU64(&output, summary.referenced_variant_ct, temporary);
  WriteU64(&output, summary.output_bytes, temporary);
  WriteU64(&output, 0, temporary);
  for (const uint64_t word : merged) WriteU64(&output, word, temporary);
  output.close();
  if (!output) throw std::runtime_error("cannot finish " + temporary);
  std::filesystem::rename(temporary, output_path);
  cleanup.Release();
  return summary;
}

struct ScoreFragmentReader::Impl {
  std::string path;
  int fd = -1;
  void* mapping = nullptr;
  uint64_t bytes = 0;
  uint64_t variant_ct = 0;
  uint32_t tile_size = 0;
  uint32_t tile_ct = 0;
  uint64_t signature_lo = 0;
  uint64_t signature_hi = 0;
  uint64_t weight_ct = 0;
  uint64_t tile_directory_offset = 0;
  uint64_t tile_data_offset = 0;
  std::vector<FragmentScore> scores;

  ~Impl() {
    if (mapping) munmap(mapping, static_cast<size_t>(bytes));
    if (fd >= 0) close(fd);
  }
};

ScoreFragmentReader::ScoreFragmentReader(const std::string& path)
    : impl_(std::make_unique<Impl>()) {
  impl_->path = path;
  impl_->fd = open(path.c_str(), O_RDONLY);
  if (impl_->fd < 0) {
    throw std::runtime_error("cannot open score fragment " + path + ": " +
                             std::strerror(errno));
  }
  impl_->bytes = ReadFileSize(impl_->fd, path);
  if (impl_->bytes < kHeaderBytes) {
    throw std::runtime_error(path + " is not a complete score fragment");
  }
  impl_->mapping = mmap(nullptr, static_cast<size_t>(impl_->bytes), PROT_READ,
                        MAP_SHARED, impl_->fd, 0);
  if (impl_->mapping == MAP_FAILED) {
    impl_->mapping = nullptr;
    throw std::runtime_error("cannot map score fragment " + path + ": " +
                             std::strerror(errno));
  }
  const auto* bytes = static_cast<const unsigned char*>(impl_->mapping);
  if (std::memcmp(bytes, kMagic, sizeof(kMagic)) != 0 ||
      GetU32(bytes + 8) != kVersion || GetU32(bytes + 12) != kHeaderBytes) {
    throw std::runtime_error(
        path + " is not a score-major score fragment; rebuild it");
  }
  impl_->variant_ct = GetU64(bytes + 16);
  impl_->tile_size = GetU32(bytes + 24);
  impl_->tile_ct = GetU32(bytes + 28);
  impl_->signature_lo = GetU64(bytes + 32);
  impl_->signature_hi = GetU64(bytes + 40);
  const uint32_t score_ct = GetU32(bytes + 48);
  impl_->weight_ct = GetU64(bytes + 56);
  const uint64_t score_records_offset = GetU64(bytes + 64);
  impl_->tile_directory_offset = GetU64(bytes + 72);
  impl_->tile_data_offset = GetU64(bytes + 80);
  const uint64_t stated_bytes = GetU64(bytes + 88);
  const uint64_t expected_tile_ct =
      impl_->tile_size
          ? (impl_->variant_ct + impl_->tile_size - 1) / impl_->tile_size
          : 0;
  const bool valid =
      impl_->variant_ct && impl_->variant_ct <= UINT32_MAX &&
      impl_->tile_size && impl_->tile_size <= kIndexMask &&
      impl_->tile_ct == expected_tile_ct && score_ct &&
      score_records_offset == kHeaderBytes &&
      score_records_offset <= impl_->tile_directory_offset &&
      impl_->tile_directory_offset <= impl_->tile_data_offset &&
      impl_->tile_data_offset <= stated_bytes && stated_bytes == impl_->bytes &&
      (static_cast<uint64_t>(impl_->tile_ct) + 1) * 8 <=
          impl_->tile_data_offset - impl_->tile_directory_offset;
  if (!valid) {
    throw std::runtime_error(path + " has invalid score-fragment dimensions");
  }
  const unsigned char* cursor = bytes + score_records_offset;
  const unsigned char* score_end = bytes + impl_->tile_directory_offset;
  impl_->scores.reserve(score_ct);
  std::unordered_set<std::string> seen_ids;
  uint64_t metadata_supported_weight_ct = 0;
  for (uint32_t idx = 0; idx < score_ct; ++idx) {
    FragmentScore score;
    cursor = ReadString(cursor, score_end, path, &score.score_id);
    if (score.score_id.empty() || !seen_ids.insert(score.score_id).second) {
      throw std::runtime_error(path + " has an invalid stable score ID");
    }
    cursor = ReadString(cursor, score_end, path, &score.info.id);
    cursor = ReadString(cursor, score_end, path, &score.info.path);
    if (score_end - cursor < 128) {
      throw std::runtime_error(path + " is truncated");
    }
    score.info.input_weight_ct = GetU64(cursor);
    score.info.zero_weight_ct = GetU64(cursor + 8);
    score.info.excluded_weight_ct = GetU64(cursor + 16);
    score.info.duplicate_weight_ct = GetU64(cursor + 24);
    score.info.catalog_weight_ct = GetU64(cursor + 32);
    score.info.matched_weight_ct = GetU64(cursor + 40);
    score.info.missing_variant_ct = GetU64(cursor + 48);
    score.info.missing_frequency_ct = GetU64(cursor + 56);
    score.info.alt_effect_ct = GetU64(cursor + 64);
    score.info.ref_effect_ct = GetU64(cursor + 72);
    const uint64_t mass_bits[] = {
        GetU64(cursor + 80), GetU64(cursor + 88), GetU64(cursor + 96),
        GetU64(cursor + 104), GetU64(cursor + 112), GetU64(cursor + 120)};
    double* masses[] = {
        &score.info.nonzero_weight_l1, &score.info.nonzero_weight_l2,
        &score.info.catalog_weight_l1, &score.info.catalog_weight_l2,
        &score.info.supported_weight_l1, &score.info.supported_weight_l2};
    for (uint32_t mass_idx = 0; mass_idx < 6; ++mass_idx) {
      std::memcpy(masses[mass_idx], &mass_bits[mass_idx], sizeof(double));
    }
    cursor += 128;
    if (score.info.catalog_weight_ct < score.info.missing_variant_ct ||
        score.info.catalog_weight_ct - score.info.missing_variant_ct <
            score.info.missing_frequency_ct ||
        score.info.matched_weight_ct != score.info.missing_frequency_ct ||
        score.info.alt_effect_ct + score.info.ref_effect_ct !=
            score.info.missing_frequency_ct) {
      throw std::runtime_error(path + " has invalid projected score counts");
    }
    metadata_supported_weight_ct +=
        score.info.catalog_weight_ct - score.info.missing_variant_ct -
        score.info.missing_frequency_ct;
    impl_->scores.push_back(std::move(score));
  }
  if (cursor != score_end) {
    throw std::runtime_error(path + " has trailing score metadata");
  }
  if (metadata_supported_weight_ct != impl_->weight_ct) {
    throw std::runtime_error(path +
                             " has inconsistent supported-weight counts");
  }
  const unsigned char* directory = bytes + impl_->tile_directory_offset;
  uint64_t previous = impl_->tile_data_offset;
  for (uint32_t idx = 0; idx <= impl_->tile_ct; ++idx) {
    const uint64_t offset = GetU64(directory + static_cast<uint64_t>(idx) * 8);
    if (offset < previous || offset > impl_->bytes ||
        (!idx && offset != impl_->tile_data_offset) ||
        (idx == impl_->tile_ct && offset != impl_->bytes)) {
      throw std::runtime_error(path + " has an invalid tile directory");
    }
    previous = offset;
  }
}

ScoreFragmentReader::~ScoreFragmentReader() = default;
ScoreFragmentReader::ScoreFragmentReader(ScoreFragmentReader&&) noexcept =
    default;
ScoreFragmentReader& ScoreFragmentReader::operator=(
    ScoreFragmentReader&&) noexcept = default;

uint64_t ScoreFragmentReader::variant_ct() const { return impl_->variant_ct; }
uint32_t ScoreFragmentReader::tile_size() const { return impl_->tile_size; }
uint32_t ScoreFragmentReader::tile_ct() const { return impl_->tile_ct; }
uint64_t ScoreFragmentReader::signature_lo() const {
  return impl_->signature_lo;
}
uint64_t ScoreFragmentReader::signature_hi() const {
  return impl_->signature_hi;
}
uint64_t ScoreFragmentReader::weight_ct() const { return impl_->weight_ct; }
uint64_t ScoreFragmentReader::file_bytes() const { return impl_->bytes; }
const std::vector<FragmentScore>& ScoreFragmentReader::scores() const {
  return impl_->scores;
}

ScoreMajorFragmentEdge ScoreFragmentScoreRow::edge(uint32_t edge_idx) const {
  if (edge_idx >= edge_ct_) {
    throw std::out_of_range("score-fragment score-row edge");
  }
  const unsigned char* cursor = edge_data_ + static_cast<uint64_t>(edge_idx) * 12;
  const uint32_t variant_and_effect = GetU32(cursor);
  const uint32_t local_variant_idx = variant_and_effect & kIndexMask;
  const uint64_t bits = GetU64(cursor + 4);
  double weight = 0.0;
  std::memcpy(&weight, &bits, sizeof(weight));
  if (local_variant_idx >= tile_variant_ct_ || !std::isfinite(weight)) {
    throw std::runtime_error(*path_ + " has an invalid score-major edge");
  }
  const bool ref_effect = variant_and_effect & kRefEffectMask;
  return {local_variant_idx, ref_effect ? -weight : weight, ref_effect};
}

void ScoreFragmentTile::OrReferencedVariants(
    std::vector<uint64_t>* words) const {
  if (words->size() != bitmap_word_ct_) {
    throw std::runtime_error("score-fragment bitmap shape mismatch");
  }
  for (uint32_t idx = 0; idx < bitmap_word_ct_; ++idx) {
    (*words)[idx] |= GetU64(bitmap_data_ + static_cast<uint64_t>(idx) * 8);
  }
}

ScoreFragmentTile ScoreFragmentReader::OpenTile(uint32_t tile_idx) const {
  if (tile_idx >= impl_->tile_ct) {
    throw std::out_of_range("score-fragment tile");
  }
  const auto* bytes = static_cast<const unsigned char*>(impl_->mapping);
  const auto* directory = bytes + impl_->tile_directory_offset;
  const uint64_t begin =
      GetU64(directory + static_cast<uint64_t>(tile_idx) * 8);
  const uint64_t end =
      GetU64(directory + static_cast<uint64_t>(tile_idx + 1) * 8);
  const unsigned char* cursor = bytes + begin;
  const unsigned char* tile_end = bytes + end;
  if (tile_end - cursor < 16) {
    throw std::runtime_error(impl_->path + " tile is truncated");
  }
  ScoreFragmentTile result;
  result.tile_idx_ = tile_idx;
  result.first_ordinal_ = tile_idx * impl_->tile_size;
  result.variant_ct_ = GetU32(cursor);
  result.referenced_variant_ct_ = GetU32(cursor + 4);
  const uint32_t score_row_ct = GetU32(cursor + 8);
  result.bitmap_word_ct_ = GetU32(cursor + 12);
  cursor += 16;
  const uint32_t expected_variant_ct = static_cast<uint32_t>(
      std::min<uint64_t>(impl_->tile_size,
                         impl_->variant_ct - result.first_ordinal_));
  const uint32_t expected_bitmap_word_ct = (expected_variant_ct + 63) / 64;
  if (result.variant_ct_ != expected_variant_ct ||
      result.bitmap_word_ct_ != expected_bitmap_word_ct ||
      static_cast<uint64_t>(tile_end - cursor) <
          static_cast<uint64_t>(result.bitmap_word_ct_) * 8) {
    throw std::runtime_error(impl_->path + " has an invalid tile header");
  }
  result.bitmap_data_ = cursor;
  uint32_t observed_referenced_ct = 0;
  for (uint32_t word_idx = 0; word_idx < result.bitmap_word_ct_; ++word_idx) {
    uint64_t word = GetU64(cursor + static_cast<uint64_t>(word_idx) * 8);
    if (word_idx + 1 == result.bitmap_word_ct_ &&
        result.variant_ct_ % 64) {
      const uint64_t valid_mask =
          (uint64_t{1} << (result.variant_ct_ % 64)) - 1;
      if (word & ~valid_mask) {
        throw std::runtime_error(impl_->path + " has bits past the tile end");
      }
    }
    observed_referenced_ct += Popcount(word);
  }
  if (observed_referenced_ct != result.referenced_variant_ct_) {
    throw std::runtime_error(impl_->path + " tile bitmap count disagrees");
  }
  cursor += static_cast<uint64_t>(result.bitmap_word_ct_) * 8;
  result.rows_.reserve(score_row_ct);
  uint32_t previous_score_idx = UINT32_MAX;
  for (uint32_t row_idx = 0; row_idx < score_row_ct; ++row_idx) {
    if (tile_end - cursor < 8) {
      throw std::runtime_error(impl_->path + " tile row is truncated");
    }
    ScoreFragmentScoreRow row;
    row.local_score_idx_ = GetU32(cursor);
    row.edge_ct_ = GetU32(cursor + 4);
    row.tile_variant_ct_ = result.variant_ct_;
    row.edge_data_ = cursor + 8;
    row.path_ = &impl_->path;
    if (row.local_score_idx_ >= impl_->scores.size() || !row.edge_ct_ ||
        (row_idx && row.local_score_idx_ <= previous_score_idx) ||
        static_cast<uint64_t>(tile_end - row.edge_data_) <
            static_cast<uint64_t>(row.edge_ct_) * 12) {
      throw std::runtime_error(impl_->path + " has an invalid tile score row");
    }
    previous_score_idx = row.local_score_idx_;
    result.rows_.push_back(row);
    cursor = row.edge_data_ + static_cast<uint64_t>(row.edge_ct_) * 12;
  }
  if (cursor != tile_end) {
    throw std::runtime_error(impl_->path + " tile has trailing data");
  }
  return result;
}

ScoreFragmentSupportSummary MeasureScoreFragmentSupport(
    const ScoreFragmentReader& fragment, const SupportIndex& support) {
  if (fragment.variant_ct() != support.variant_ct() ||
      fragment.signature_lo() != support.signature_lo() ||
      fragment.signature_hi() != support.signature_hi()) {
    throw std::runtime_error(
        "score fragment and support index describe different variant indexes");
  }
  ScoreFragmentSupportSummary result;
  result.variant_ct = fragment.variant_ct();
  result.score_ct = fragment.scores().size();
  result.scores.reserve(fragment.scores().size());
  for (const auto& score : fragment.scores()) {
    ScoreFragmentSupportRow row;
    row.score_id = score.score_id;
    row.column_name = score.info.id;
    row.reference_weight_ct =
        score.info.catalog_weight_ct - score.info.missing_variant_ct -
        score.info.missing_frequency_ct;
    row.reference_weight_l1 = score.info.supported_weight_l1;
    row.reference_weight_l2_squared = score.info.supported_weight_l2;
    result.reference_weight_ct += row.reference_weight_ct;
    result.scores.push_back(std::move(row));
  }
  std::vector<uint64_t> observed_reference_counts(result.scores.size(), 0);
  std::vector<double> observed_reference_l1(result.scores.size(), 0.0);
  std::vector<double> observed_reference_l2_squared(result.scores.size(), 0.0);
  for (uint32_t tile_idx = 0; tile_idx < fragment.tile_ct(); ++tile_idx) {
    const auto tile = fragment.OpenTile(tile_idx);
    for (const auto& score_row : tile.rows()) {
      auto& output = result.scores.at(score_row.local_score_idx());
      for (uint32_t edge_idx = 0; edge_idx < score_row.edge_ct(); ++edge_idx) {
        const auto edge = score_row.edge(edge_idx);
        const uint32_t ordinal = tile.first_ordinal() + edge.local_variant_idx;
        const double weight = std::abs(edge.beta_alt);
        ++observed_reference_counts[score_row.local_score_idx()];
        observed_reference_l1[score_row.local_score_idx()] += weight;
        observed_reference_l2_squared[score_row.local_score_idx()] +=
            weight * weight;
        if (support.state(ordinal) != VariantSupport::kUsable) continue;
        ++output.available_weight_ct;
        output.available_weight_l1 += weight;
        output.available_weight_l2_squared += weight * weight;
        ++result.available_weight_ct;
      }
    }
  }
  for (uint32_t score_idx = 0; score_idx < result.scores.size(); ++score_idx) {
    const auto& row = result.scores[score_idx];
    auto close = [](double observed, double expected) {
      return std::abs(observed - expected) <=
             1e-10 * std::max(1.0, std::abs(expected));
    };
    if (row.available_weight_ct > row.reference_weight_ct ||
        observed_reference_counts[score_idx] != row.reference_weight_ct ||
        !close(observed_reference_l1[score_idx], row.reference_weight_l1) ||
        !close(observed_reference_l2_squared[score_idx],
               row.reference_weight_l2_squared)) {
      throw std::runtime_error(
          "score-fragment support accounting is inconsistent for " +
          row.score_id);
    }
  }
  if (result.reference_weight_ct != fragment.weight_ct()) {
    throw std::runtime_error(
        "score-fragment support report has inconsistent reference weights");
  }
  return result;
}

}  // namespace pgensparsescore
