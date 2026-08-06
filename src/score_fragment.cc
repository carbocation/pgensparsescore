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
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "io.h"

namespace pgensparsescore {
namespace {

constexpr char kMagic[8] = {'P', 'G', 'S', 'S', 'F', 'R', 'A', 'G'};
constexpr uint32_t kVersion = 1;
constexpr uint32_t kHeaderBytes = 96;
constexpr uint32_t kRefEffectMask = 0x80000000U;
constexpr uint32_t kScoreIndexMask = 0x7fffffffU;
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
  const auto timestamp = std::chrono::steady_clock::now()
                             .time_since_epoch()
                             .count();
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
    result.push_back({fields[score_id_idx], fields[column_idx], fields[path_idx],
                      resolved.lexically_normal().string()});
  }
  if (result.empty()) throw std::runtime_error(path + " has no scores");
  return result;
}

void FillHeader(unsigned char* header, const VariantIndex& index,
                uint32_t score_ct, uint64_t weight_ct,
                uint64_t score_records_offset,
                uint64_t block_directory_offset, uint64_t block_data_offset,
                uint64_t file_bytes) {
  std::memset(header, 0, kHeaderBytes);
  std::memcpy(header, kMagic, sizeof(kMagic));
  PutU32(header + 8, kVersion);
  PutU32(header + 12, kHeaderBytes);
  PutU64(header + 16, index.variant_ct());
  PutU32(header + 24, index.block_size());
  PutU32(header + 28, index.block_ct());
  PutU64(header + 32, index.signature_lo());
  PutU64(header + 40, index.signature_hi());
  PutU32(header + 48, score_ct);
  PutU64(header + 56, weight_ct);
  PutU64(header + 64, score_records_offset);
  PutU64(header + 72, block_directory_offset);
  PutU64(header + 80, block_data_offset);
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

}  // namespace

ScoreFragmentSummary CompileScoreFragment(
    const ScoreFragmentCompileOptions& options, ProgressReporter* progress) {
  if (std::filesystem::exists(options.output_path)) {
    throw std::runtime_error("score-fragment output already exists: " +
                             options.output_path);
  }
  VariantIndex index(options.variant_index_path);
  const auto manifest = ReadManifest(options.manifest_path);
  if (manifest.size() > kScoreIndexMask) {
    throw std::runtime_error("score fragment contains too many scores");
  }
  if (progress) {
    progress->Event("fragment", "start",
                    {{"score_files_total", manifest.size()},
                     {"variant_index_variants", index.variant_ct()},
                     {"variant_blocks", index.block_ct()}},
                    {{"manifest", options.manifest_path},
                     {"variant_index", options.variant_index_path},
                     {"output", options.output_path}});
  }

  const auto temp_directory = CreateTemporaryDirectory(options);
  DirectoryCleanup cleanup(temp_directory);
  std::vector<std::unique_ptr<std::ofstream>> buckets(index.block_ct());
  std::vector<uint64_t> bucket_weight_ct(index.block_ct(), 0);
  std::vector<FragmentScore> scores;
  scores.reserve(manifest.size());
  ScoreFragmentSummary summary;
  summary.variant_index_variant_ct = index.variant_ct();
  summary.block_size = index.block_size();
  summary.block_ct = index.block_ct();
  summary.score_ct = static_cast<uint32_t>(manifest.size());

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
      BucketWeight record{*ordinal,
                          manifest_idx |
                              (ref_effect ? kRefEffectMask : 0U),
                          weight};
      buckets[block_idx]->write(reinterpret_cast<const char*>(&record),
                                sizeof(record));
      if (!*buckets[block_idx]) {
        throw std::runtime_error("cannot write score-fragment bucket");
      }
      ++bucket_weight_ct[block_idx];
      ++info.catalog_weight_ct;
      ++summary.weight_ct;
      if (progress && !(summary.input_weight_ct % 1000000)) {
        progress->MaybeEvent(
            "fragment", "parse_weights",
            {{"score_files_processed", manifest_idx},
             {"score_files_total", manifest.size()},
             {"input_weights", summary.input_weight_ct},
             {"retained_weights", summary.weight_ct},
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
           {"retained_weights", summary.weight_ct},
           {"zero_weights", summary.zero_weight_ct},
           {"excluded_weights", summary.excluded_weight_ct},
           {"duplicate_weights", summary.duplicate_weight_ct}},
          {{"current_score", item.score_id}});
    }
  }
  for (auto& bucket : buckets) {
    if (bucket) {
      bucket->close();
      if (!*bucket) throw std::runtime_error("cannot finish fragment bucket");
      bucket.reset();
    }
  }

  const std::string temporary_output = options.output_path + ".tmp";
  if (std::filesystem::exists(temporary_output)) {
    throw std::runtime_error("temporary fragment output already exists: " +
                             temporary_output);
  }
  FileCleanup output_cleanup(temporary_output);
  std::ofstream output(temporary_output,
                       std::ios::binary | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot create " + temporary_output);
  std::vector<unsigned char> header(kHeaderBytes, 0);
  WriteBytes(&output, header.data(), header.size(), temporary_output);
  const uint64_t score_records_offset = Position(&output, temporary_output);
  for (const auto& score : scores) {
    WriteString(&output, score.score_id, temporary_output);
    WriteString(&output, score.info.id, temporary_output);
    WriteString(&output, score.info.path, temporary_output);
    WriteU64(&output, score.info.input_weight_ct, temporary_output);
    WriteU64(&output, score.info.zero_weight_ct, temporary_output);
    WriteU64(&output, score.info.excluded_weight_ct, temporary_output);
    WriteU64(&output, score.info.duplicate_weight_ct, temporary_output);
    WriteU64(&output, score.info.catalog_weight_ct, temporary_output);
  }
  const uint64_t block_directory_offset = Position(&output, temporary_output);
  std::vector<uint64_t> block_offsets(index.block_ct() + 1, 0);
  for (uint32_t idx = 0; idx <= index.block_ct(); ++idx) {
    WriteU64(&output, 0, temporary_output);
  }
  const uint64_t block_data_offset = Position(&output, temporary_output);
  uint64_t weights_written = 0;
  for (uint32_t block_idx = 0; block_idx < index.block_ct(); ++block_idx) {
    block_offsets[block_idx] = Position(&output, temporary_output);
    std::vector<BucketWeight> records(bucket_weight_ct[block_idx]);
    if (!records.empty()) {
      const auto path = temp_directory /
                        ("block-" + std::to_string(block_idx) + ".bin");
      std::ifstream input(path, std::ios::binary);
      if (!input) throw std::runtime_error("cannot open " + path.string());
      input.read(reinterpret_cast<char*>(records.data()),
                 static_cast<std::streamsize>(records.size() * sizeof(records[0])));
      if (!input || input.peek() != std::ifstream::traits_type::eof()) {
        throw std::runtime_error("invalid score-fragment bucket " + path.string());
      }
      std::sort(records.begin(), records.end(),
                [](const BucketWeight& lhs, const BucketWeight& rhs) {
                  if (lhs.ordinal != rhs.ordinal) return lhs.ordinal < rhs.ordinal;
                  return lhs.score_and_effect < rhs.score_and_effect;
                });
    }
    uint32_t variant_group_ct = 0;
    for (size_t idx = 0; idx < records.size(); ++idx) {
      if (!idx || records[idx].ordinal != records[idx - 1].ordinal) {
        ++variant_group_ct;
      }
    }
    WriteU32(&output, variant_group_ct, temporary_output);
    size_t begin = 0;
    while (begin < records.size()) {
      size_t end = begin + 1;
      while (end < records.size() &&
             records[end].ordinal == records[begin].ordinal) {
        ++end;
      }
      const uint32_t ordinal = records[begin].ordinal;
      const uint32_t expected_block = ordinal / index.block_size();
      if (expected_block != block_idx || end - begin > UINT32_MAX) {
        throw std::runtime_error("invalid score-fragment bucket contents");
      }
      WriteU32(&output, ordinal, temporary_output);
      WriteU32(&output, static_cast<uint32_t>(end - begin), temporary_output);
      for (size_t idx = begin; idx < end; ++idx) {
        WriteU32(&output, records[idx].score_and_effect, temporary_output);
        WriteDouble(&output, records[idx].weight, temporary_output);
      }
      weights_written += end - begin;
      begin = end;
    }
    if (progress) {
      progress->MaybeEvent("fragment", "serialize_blocks",
                           {{"blocks_written", block_idx + 1},
                            {"blocks_total", index.block_ct()},
                            {"weights_written", weights_written},
                            {"weights_total", summary.weight_ct},
                            {"output_bytes", Position(&output, temporary_output)}});
    }
  }
  block_offsets[index.block_ct()] = Position(&output, temporary_output);
  if (weights_written != summary.weight_ct) {
    throw std::runtime_error("fragment weight count changed during serialization");
  }
  summary.output_bytes = block_offsets.back();
  output.seekp(static_cast<std::streamoff>(block_directory_offset));
  for (const uint64_t offset : block_offsets) {
    WriteU64(&output, offset, temporary_output);
  }
  FillHeader(header.data(), index, scores.size(), summary.weight_ct,
             score_records_offset, block_directory_offset, block_data_offset,
             summary.output_bytes);
  output.seekp(0);
  WriteBytes(&output, header.data(), header.size(), temporary_output);
  output.close();
  if (!output) throw std::runtime_error("cannot finish " + temporary_output);
  std::filesystem::rename(temporary_output, options.output_path);
  output_cleanup.Release();
  if (progress) {
    progress->Event("fragment", "complete",
                    {{"scores", summary.score_ct},
                     {"input_weights", summary.input_weight_ct},
                     {"retained_weights", summary.weight_ct},
                     {"zero_weights", summary.zero_weight_ct},
                     {"excluded_weights", summary.excluded_weight_ct},
                     {"duplicate_weights", summary.duplicate_weight_ct},
                     {"output_bytes", summary.output_bytes}});
  }
  return summary;
}

struct ScoreFragmentReader::Impl {
  std::string path;
  int fd = -1;
  void* mapping = nullptr;
  uint64_t bytes = 0;
  uint64_t variant_ct = 0;
  uint32_t block_size = 0;
  uint32_t block_ct = 0;
  uint64_t signature_lo = 0;
  uint64_t signature_hi = 0;
  uint64_t weight_ct = 0;
  uint64_t block_directory_offset = 0;
  uint64_t block_data_offset = 0;
  std::vector<FragmentScore> scores;
  std::unordered_set<uint32_t> score_indexes;

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
    throw std::runtime_error(path + " has an invalid score-fragment header");
  }
  impl_->variant_ct = GetU64(bytes + 16);
  impl_->block_size = GetU32(bytes + 24);
  impl_->block_ct = GetU32(bytes + 28);
  impl_->signature_lo = GetU64(bytes + 32);
  impl_->signature_hi = GetU64(bytes + 40);
  const uint32_t score_ct = GetU32(bytes + 48);
  impl_->weight_ct = GetU64(bytes + 56);
  const uint64_t score_records_offset = GetU64(bytes + 64);
  impl_->block_directory_offset = GetU64(bytes + 72);
  impl_->block_data_offset = GetU64(bytes + 80);
  const uint64_t stated_bytes = GetU64(bytes + 88);
  const uint64_t expected_block_ct =
      impl_->block_size
          ? (impl_->variant_ct + impl_->block_size - 1) / impl_->block_size
          : 0;
  const bool valid = impl_->variant_ct && impl_->variant_ct <= UINT32_MAX &&
                     impl_->block_size && impl_->block_ct == expected_block_ct &&
                     score_ct && score_records_offset == kHeaderBytes &&
                     score_records_offset <= impl_->block_directory_offset &&
                     impl_->block_directory_offset <= impl_->block_data_offset &&
                     impl_->block_data_offset <= stated_bytes &&
                     stated_bytes == impl_->bytes &&
                     static_cast<uint64_t>(impl_->block_ct + 1) * 8 <=
                         impl_->block_data_offset - impl_->block_directory_offset;
  if (!valid) {
    throw std::runtime_error(path + " has invalid score-fragment dimensions");
  }
  const unsigned char* cursor = bytes + score_records_offset;
  const unsigned char* score_end = bytes + impl_->block_directory_offset;
  impl_->scores.reserve(score_ct);
  std::unordered_set<std::string> seen_ids;
  for (uint32_t idx = 0; idx < score_ct; ++idx) {
    FragmentScore score;
    cursor = ReadString(cursor, score_end, path, &score.score_id);
    if (score.score_id.empty() || !seen_ids.insert(score.score_id).second) {
      throw std::runtime_error(path + " has an invalid stable score ID");
    }
    impl_->score_indexes.insert(idx);
    cursor = ReadString(cursor, score_end, path, &score.info.id);
    cursor = ReadString(cursor, score_end, path, &score.info.path);
    if (score_end - cursor < 40) throw std::runtime_error(path + " is truncated");
    score.info.input_weight_ct = GetU64(cursor);
    score.info.zero_weight_ct = GetU64(cursor + 8);
    score.info.excluded_weight_ct = GetU64(cursor + 16);
    score.info.duplicate_weight_ct = GetU64(cursor + 24);
    score.info.catalog_weight_ct = GetU64(cursor + 32);
    cursor += 40;
    impl_->scores.push_back(std::move(score));
  }
  if (cursor != score_end) {
    throw std::runtime_error(path + " has trailing score metadata");
  }
  const unsigned char* directory = bytes + impl_->block_directory_offset;
  uint64_t previous = impl_->block_data_offset;
  for (uint32_t idx = 0; idx <= impl_->block_ct; ++idx) {
    const uint64_t offset = GetU64(directory + static_cast<uint64_t>(idx) * 8);
    if (offset < previous || offset > impl_->bytes ||
        (!idx && offset != impl_->block_data_offset) ||
        (idx == impl_->block_ct && offset != impl_->bytes)) {
      throw std::runtime_error(path + " has an invalid block directory");
    }
    previous = offset;
  }
}

ScoreFragmentReader::~ScoreFragmentReader() = default;
ScoreFragmentReader::ScoreFragmentReader(ScoreFragmentReader&&) noexcept = default;
ScoreFragmentReader& ScoreFragmentReader::operator=(ScoreFragmentReader&&) noexcept =
    default;

uint64_t ScoreFragmentReader::variant_ct() const { return impl_->variant_ct; }
uint32_t ScoreFragmentReader::block_size() const { return impl_->block_size; }
uint32_t ScoreFragmentReader::block_ct() const { return impl_->block_ct; }
uint64_t ScoreFragmentReader::signature_lo() const { return impl_->signature_lo; }
uint64_t ScoreFragmentReader::signature_hi() const { return impl_->signature_hi; }
uint64_t ScoreFragmentReader::weight_ct() const { return impl_->weight_ct; }
uint64_t ScoreFragmentReader::file_bytes() const { return impl_->bytes; }
const std::vector<FragmentScore>& ScoreFragmentReader::scores() const {
  return impl_->scores;
}

struct ScoreFragmentBlockCursor::Impl {
  const unsigned char* cursor = nullptr;
  const unsigned char* end = nullptr;
  const unsigned char* edge_data = nullptr;
  const std::unordered_set<uint32_t>* score_indexes = nullptr;
  const std::string* path = nullptr;
  uint32_t groups_remaining = 0;
  uint32_t current_ordinal = 0;
  uint32_t previous_ordinal = 0;
  uint32_t edge_ct = 0;
  uint32_t block_idx = 0;
  uint32_t block_size = 0;
  uint64_t variant_ct = 0;
  bool has_previous = false;

  void ParseCurrent() {
    if (!groups_remaining) {
      if (cursor != end) {
        throw std::runtime_error(*path + " block has trailing data");
      }
      return;
    }
    if (end - cursor < 8) {
      throw std::runtime_error(*path + " block is truncated");
    }
    current_ordinal = GetU32(cursor);
    edge_ct = GetU32(cursor + 4);
    edge_data = cursor + 8;
    if (!edge_ct || current_ordinal >= variant_ct ||
        current_ordinal / block_size != block_idx ||
        (has_previous && current_ordinal <= previous_ordinal) ||
        static_cast<uint64_t>(end - edge_data) <
            static_cast<uint64_t>(edge_ct) * 12) {
      throw std::runtime_error(*path + " has an invalid block record");
    }
  }
};

ScoreFragmentBlockCursor::ScoreFragmentBlockCursor()
    : impl_(std::make_unique<Impl>()) {}
ScoreFragmentBlockCursor::~ScoreFragmentBlockCursor() = default;
ScoreFragmentBlockCursor::ScoreFragmentBlockCursor(
    ScoreFragmentBlockCursor&&) noexcept = default;
ScoreFragmentBlockCursor& ScoreFragmentBlockCursor::operator=(
    ScoreFragmentBlockCursor&&) noexcept = default;

bool ScoreFragmentBlockCursor::done() const {
  return !impl_->groups_remaining;
}

uint32_t ScoreFragmentBlockCursor::ordinal() const {
  if (done()) throw std::out_of_range("score-fragment block cursor");
  return impl_->current_ordinal;
}

void ScoreFragmentBlockCursor::AppendEdges(std::vector<Edge>* edges) const {
  if (done()) throw std::out_of_range("score-fragment block cursor");
  const unsigned char* cursor = impl_->edge_data;
  edges->reserve(edges->size() + impl_->edge_ct);
  for (uint32_t edge_idx = 0; edge_idx < impl_->edge_ct; ++edge_idx) {
    const uint32_t score_and_effect = GetU32(cursor);
    const uint32_t score_idx = score_and_effect & kScoreIndexMask;
    const uint64_t bits = GetU64(cursor + 4);
    double weight = 0.0;
    std::memcpy(&weight, &bits, sizeof(weight));
    if (!std::isfinite(weight) || !impl_->score_indexes->count(score_idx)) {
      throw std::runtime_error(*impl_->path +
                               " has an invalid fragment weight record");
    }
    const bool ref_effect = score_and_effect & kRefEffectMask;
    edges->push_back(
        {score_idx, ref_effect ? -weight : weight, ref_effect});
    cursor += 12;
  }
}

void ScoreFragmentBlockCursor::Next() {
  if (done()) throw std::out_of_range("score-fragment block cursor");
  impl_->cursor = impl_->edge_data + static_cast<uint64_t>(impl_->edge_ct) * 12;
  impl_->previous_ordinal = impl_->current_ordinal;
  impl_->has_previous = true;
  --impl_->groups_remaining;
  impl_->ParseCurrent();
}

ScoreFragmentBlockCursor ScoreFragmentReader::OpenBlock(
    uint32_t block_idx) const {
  if (block_idx >= impl_->block_ct) throw std::out_of_range("fragment block");
  const auto* bytes = static_cast<const unsigned char*>(impl_->mapping);
  const auto* directory = bytes + impl_->block_directory_offset;
  const uint64_t begin = GetU64(directory + static_cast<uint64_t>(block_idx) * 8);
  const uint64_t end =
      GetU64(directory + static_cast<uint64_t>(block_idx + 1) * 8);
  ScoreFragmentBlockCursor result;
  result.impl_->cursor = bytes + begin;
  result.impl_->end = bytes + end;
  result.impl_->score_indexes = &impl_->score_indexes;
  result.impl_->path = &impl_->path;
  result.impl_->block_idx = block_idx;
  result.impl_->block_size = impl_->block_size;
  result.impl_->variant_ct = impl_->variant_ct;
  if (result.impl_->end - result.impl_->cursor < 4) {
    throw std::runtime_error(impl_->path + " block is truncated");
  }
  result.impl_->groups_remaining = GetU32(result.impl_->cursor);
  result.impl_->cursor += 4;
  result.impl_->ParseCurrent();
  return result;
}

void ScoreFragmentReader::ReadBlock(
    uint32_t block_idx, std::vector<IndexedVariantEdges>* variants) const {
  variants->clear();
  auto cursor = OpenBlock(block_idx);
  while (!cursor.done()) {
    IndexedVariantEdges variant;
    variant.ordinal = cursor.ordinal();
    cursor.AppendEdges(&variant.edges);
    variants->push_back(std::move(variant));
    cursor.Next();
  }
}

}  // namespace pgensparsescore
