// SPDX-License-Identifier: GPL-3.0-only
#include "compiled_catalog.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "io.h"

namespace pgensparsescore {
namespace {

using Header = std::unordered_map<std::string, size_t>;
constexpr char kMagic[] = {'P', 'G', 'S', 'S', 'C', 'A', 'T', '1'};
constexpr uint32_t kFormatVersion = 2;
constexpr uint32_t kMaximumStringBytes = 1U << 28;

Header MakeHeader(const std::vector<std::string>& fields) {
  Header result;
  for (size_t idx = 0; idx < fields.size(); ++idx) {
    std::string name = fields[idx];
    if (!name.empty() && name.front() == '#') {
      name.erase(0, 1);
    }
    if (!result.emplace(name, idx).second) {
      throw std::runtime_error("duplicate column in header: " + name);
    }
  }
  return result;
}

size_t RequireColumn(const Header& header, const std::string& name,
                     const std::string& path) {
  const auto iter = header.find(name);
  if (iter == header.end()) {
    throw std::runtime_error(path + " is missing required column " + name);
  }
  return iter->second;
}

size_t FindScoreColumn(const Header& header, const std::string& path) {
  for (const char* name : {"COLUMN_NAME", "SCORE", "SCORE_ID"}) {
    const auto iter = header.find(name);
    if (iter != header.end()) {
      return iter->second;
    }
  }
  throw std::runtime_error(path +
                           " is missing COLUMN_NAME, SCORE, or SCORE_ID");
}

void RequireFields(const std::vector<std::string>& fields, size_t index,
                   const std::string& path, uint64_t line_number) {
  if (fields.size() <= index) {
    throw std::runtime_error(path + ": line " +
                             std::to_string(line_number) +
                             " has too few fields");
  }
}

double ParseWeight(const std::string& value, const std::string& path,
                   uint64_t line_number) {
  errno = 0;
  char* end = nullptr;
  const double result = std::strtod(value.c_str(), &end);
  if (errno || end == value.c_str() || *end != '\0' ||
      !std::isfinite(result)) {
    throw std::runtime_error(path + ": line " +
                             std::to_string(line_number) +
                             " has invalid finite weight: " + value);
  }
  return result;
}

void WriteBytes(std::ostream* output, const void* bytes, size_t byte_ct) {
  output->write(static_cast<const char*>(bytes),
                static_cast<std::streamsize>(byte_ct));
  if (!*output) {
    throw std::runtime_error("cannot write compiled catalog");
  }
}

void WriteU8(std::ostream* output, uint8_t value) {
  WriteBytes(output, &value, 1);
}

void WriteU32(std::ostream* output, uint32_t value) {
  uint8_t bytes[4];
  for (uint32_t idx = 0; idx < 4; ++idx) {
    bytes[idx] = static_cast<uint8_t>(value >> (8 * idx));
  }
  WriteBytes(output, bytes, sizeof(bytes));
}

void WriteU64(std::ostream* output, uint64_t value) {
  uint8_t bytes[8];
  for (uint32_t idx = 0; idx < 8; ++idx) {
    bytes[idx] = static_cast<uint8_t>(value >> (8 * idx));
  }
  WriteBytes(output, bytes, sizeof(bytes));
}

void WriteDouble(std::ostream* output, double value) {
  uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  WriteU64(output, bits);
}

void WriteString(std::ostream* output, const std::string& value) {
  if (value.size() > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error("compiled catalog string is too long");
  }
  WriteU32(output, static_cast<uint32_t>(value.size()));
  WriteBytes(output, value.data(), value.size());
}

void ReadBytes(std::istream* input, void* bytes, size_t byte_ct,
               const std::string& path) {
  input->read(static_cast<char*>(bytes), static_cast<std::streamsize>(byte_ct));
  if (!*input) {
    throw std::runtime_error(path + " is truncated");
  }
}

uint8_t ReadU8(std::istream* input, const std::string& path) {
  uint8_t value = 0;
  ReadBytes(input, &value, 1, path);
  return value;
}

uint32_t ReadU32(std::istream* input, const std::string& path) {
  uint8_t bytes[4];
  ReadBytes(input, bytes, sizeof(bytes), path);
  uint32_t value = 0;
  for (uint32_t idx = 0; idx < 4; ++idx) {
    value |= static_cast<uint32_t>(bytes[idx]) << (8 * idx);
  }
  return value;
}

uint64_t ReadU64(std::istream* input, const std::string& path) {
  uint8_t bytes[8];
  ReadBytes(input, bytes, sizeof(bytes), path);
  uint64_t value = 0;
  for (uint32_t idx = 0; idx < 8; ++idx) {
    value |= static_cast<uint64_t>(bytes[idx]) << (8 * idx);
  }
  return value;
}

double ReadDouble(std::istream* input, const std::string& path) {
  const uint64_t bits = ReadU64(input, path);
  double value = 0.0;
  std::memcpy(&value, &bits, sizeof(value));
  if (!std::isfinite(value)) {
    throw std::runtime_error(path + " contains a nonfinite weight");
  }
  return value;
}

std::string ReadString(std::istream* input, const std::string& path) {
  const uint32_t byte_ct = ReadU32(input, path);
  if (byte_ct > kMaximumStringBytes) {
    throw std::runtime_error(path + " contains an unreasonable string length");
  }
  std::string result(byte_ct, '\0');
  ReadBytes(input, result.data(), byte_ct, path);
  return result;
}

}  // namespace

CompiledCatalog CompileSourceCatalog(
    const std::string& manifest_path,
    const std::unordered_set<std::string>* included_source_ids,
    ProgressReporter* progress) {
  struct ManifestRow {
    std::string score;
    std::string display_path;
    std::string resolved_path;
  };
  std::vector<ManifestRow> manifest;
  {
    LineReader reader(manifest_path);
    std::string line;
    if (!reader.GetLine(&line)) {
      throw std::runtime_error(manifest_path + " is empty");
    }
    const Header header = MakeHeader(SplitTabs(line));
    const size_t score_idx = FindScoreColumn(header, manifest_path);
    const size_t path_idx = RequireColumn(header, "PATH", manifest_path);
    const std::filesystem::path base =
        std::filesystem::absolute(manifest_path).parent_path();
    std::unordered_set<std::string> seen_scores;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) {
        continue;
      }
      const auto fields = SplitTabs(line);
      RequireFields(fields, std::max(score_idx, path_idx), manifest_path,
                    line_number);
      if (!seen_scores.insert(fields[score_idx]).second) {
        throw std::runtime_error("duplicate manifest SCORE: " +
                                 fields[score_idx]);
      }
      std::filesystem::path resolved(fields[path_idx]);
      if (resolved.is_relative()) {
        resolved = base / resolved;
      }
      manifest.push_back({fields[score_idx], fields[path_idx],
                          resolved.lexically_normal().string()});
    }
  }
  if (manifest.empty()) {
    throw std::runtime_error(manifest_path + " has no scores");
  }
  if (progress) {
    progress->Event("compile", "manifest_loaded",
                    {{"score_files_total", manifest.size()},
                     {"included_variants",
                      included_source_ids ? included_source_ids->size() : 0}});
  }

  CompiledCatalog result;
  result.scores.reserve(manifest.size());
  std::unordered_map<std::string, uint32_t> variant_by_source;
  uint64_t input_weight_ct = 0;
  uint64_t zero_weight_ct = 0;
  uint64_t excluded_weight_ct = 0;
  uint64_t duplicate_weight_ct = 0;

  for (uint32_t score_idx = 0; score_idx < manifest.size(); ++score_idx) {
    const auto& item = manifest[score_idx];
    ScoreInfo info;
    info.id = item.score;
    info.path = item.display_path;
    LineReader reader(item.resolved_path);
    std::string line;
    if (!reader.GetLine(&line)) {
      throw std::runtime_error(item.resolved_path + " is empty");
    }
    const Header header = MakeHeader(SplitTabs(line));
    const size_t snp_idx = RequireColumn(header, "SNP", item.resolved_path);
    const size_t effect_idx =
        RequireColumn(header, "EFFECT_ALLELE", item.resolved_path);
    const size_t other_idx =
        RequireColumn(header, "OTHER_ALLELE", item.resolved_path);
    const size_t weight_idx =
        RequireColumn(header, "EFFECT_ALLELE_WEIGHT", item.resolved_path);
    const size_t max_idx =
        std::max({snp_idx, effect_idx, other_idx, weight_idx});
    std::unordered_set<std::string> score_variants;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) {
        continue;
      }
      ++info.input_weight_ct;
      ++input_weight_ct;
      if (progress && !(input_weight_ct % 1000000)) {
        progress->MaybeEvent(
            "compile", "parse_weights",
            {{"score_files_processed", score_idx},
             {"score_files_total", manifest.size()},
             {"input_weights", input_weight_ct},
             {"zero_weights", zero_weight_ct},
             {"excluded_weights", excluded_weight_ct},
             {"duplicate_weights", duplicate_weight_ct},
             {"retained_weights", result.weight_ct},
             {"unique_variants", result.variants.size()}},
            {{"current_score", item.score}});
      }
      const auto fields = SplitTabs(line);
      RequireFields(fields, max_idx, item.resolved_path, line_number);
      const double weight =
          ParseWeight(fields[weight_idx], item.resolved_path, line_number);
      if (weight == 0.0) {
        ++info.zero_weight_ct;
        ++zero_weight_ct;
        continue;
      }
      const std::string& source_id = fields[snp_idx];
      if (included_source_ids && !included_source_ids->count(source_id)) {
        ++info.excluded_weight_ct;
        ++excluded_weight_ct;
        continue;
      }
      const std::string& effect = fields[effect_idx];
      const std::string& other = fields[other_idx];
      if (source_id.empty() || effect.empty() || other.empty() ||
          effect == other) {
        throw std::runtime_error(item.resolved_path + ": line " +
                                 std::to_string(line_number) +
                                 " has invalid variant ID or allele pair");
      }
      if (!score_variants.insert(source_id).second) {
        ++info.duplicate_weight_ct;
        ++duplicate_weight_ct;
      }
      const std::string allele0 = std::min(effect, other);
      const std::string allele1 = std::max(effect, other);
      uint32_t variant_idx = 0;
      const auto existing = variant_by_source.find(source_id);
      if (existing == variant_by_source.end()) {
        if (result.variants.size() == std::numeric_limits<uint32_t>::max()) {
          throw std::runtime_error("compiled catalog exceeds variant limit");
        }
        variant_idx = static_cast<uint32_t>(result.variants.size());
        variant_by_source.emplace(source_id, variant_idx);
        result.variants.push_back({source_id, allele0, allele1, {}});
      } else {
        variant_idx = existing->second;
        const auto& variant = result.variants[variant_idx];
        if (variant.allele0 != allele0 || variant.allele1 != allele1) {
          throw std::runtime_error(
              item.resolved_path + ": line " + std::to_string(line_number) +
              " has an allele pair inconsistent with earlier rows for " +
              source_id);
        }
      }
      result.variants[variant_idx].weights.push_back(
          {score_idx, weight,
           static_cast<uint8_t>(effect == allele1 ? 1 : 0)});
      ++info.catalog_weight_ct;
      ++result.weight_ct;
    }
    result.scores.push_back(std::move(info));
    if (progress) {
      progress->MaybeEvent(
          "compile", "parse_weights",
          {{"score_files_processed", score_idx + 1},
           {"score_files_total", manifest.size()},
           {"input_weights", input_weight_ct},
           {"zero_weights", zero_weight_ct},
           {"excluded_weights", excluded_weight_ct},
           {"duplicate_weights", duplicate_weight_ct},
           {"retained_weights", result.weight_ct},
           {"unique_variants", result.variants.size()}},
          {{"current_score", item.score}});
    }
  }

  if (progress) {
    progress->Event(
        "compile", "sort_weights",
        {{"score_files_processed", result.scores.size()},
         {"score_files_total", manifest.size()},
         {"input_weights", input_weight_ct},
         {"zero_weights", zero_weight_ct},
         {"excluded_weights", excluded_weight_ct},
         {"duplicate_weights", duplicate_weight_ct},
         {"retained_weights", result.weight_ct},
         {"unique_variants", result.variants.size()}});
  }
  uint64_t variants_sorted = 0;
  uint64_t weights_sorted = 0;
  for (auto& variant : result.variants) {
    std::sort(variant.weights.begin(), variant.weights.end(),
              [](const CompiledWeight& lhs, const CompiledWeight& rhs) {
                return lhs.score_idx < rhs.score_idx;
              });
    ++variants_sorted;
    weights_sorted += variant.weights.size();
    if (progress && !(variants_sorted % 1000000)) {
      progress->MaybeEvent(
          "compile", "sort_weights",
          {{"variants_sorted", variants_sorted},
           {"variants_total", result.variants.size()},
           {"weights_sorted", weights_sorted},
           {"weights_total", result.weight_ct}});
    }
  }
  if (progress) {
    progress->Event("compile", "sort_variants",
                    {{"variants_sorted", variants_sorted},
                     {"variants_total", result.variants.size()},
                     {"weights_sorted", weights_sorted},
                     {"weights_total", result.weight_ct}});
  }
  std::sort(result.variants.begin(), result.variants.end(),
            [](const CompiledVariant& lhs, const CompiledVariant& rhs) {
              return lhs.source_id < rhs.source_id;
            });
  if (progress) {
    progress->Event("compile", "weights_ready",
                    {{"score_files_processed", result.scores.size()},
                     {"input_weights", input_weight_ct},
                     {"zero_weights", zero_weight_ct},
                     {"excluded_weights", excluded_weight_ct},
                     {"duplicate_weights", duplicate_weight_ct},
                     {"retained_weights", result.weight_ct},
                     {"unique_variants", result.variants.size()}});
  }
  return result;
}

void WriteCompiledCatalog(const std::string& path,
                          const CompiledCatalog& catalog,
                          ProgressReporter* progress) {
  const std::string temporary = path + ".tmp";
  std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
  if (!output) {
    throw std::runtime_error("cannot create " + temporary);
  }
  WriteBytes(&output, kMagic, sizeof(kMagic));
  WriteU32(&output, kFormatVersion);
  WriteU32(&output, static_cast<uint32_t>(catalog.scores.size()));
  WriteU64(&output, catalog.variants.size());
  WriteU64(&output, catalog.weight_ct);
  for (const auto& score : catalog.scores) {
    WriteString(&output, score.id);
    WriteString(&output, score.path);
    WriteU64(&output, score.input_weight_ct);
    WriteU64(&output, score.zero_weight_ct);
    WriteU64(&output, score.excluded_weight_ct);
    WriteU64(&output, score.duplicate_weight_ct);
    WriteU64(&output, score.catalog_weight_ct);
  }
  if (progress) {
    progress->Event("compile", "serialize",
                    {{"scores_written", catalog.scores.size()},
                     {"variants_written", 0},
                     {"variants_total", catalog.variants.size()},
                     {"weights_written", 0},
                     {"weights_total", catalog.weight_ct}});
  }
  uint64_t variants_written = 0;
  uint64_t weights_written = 0;
  for (const auto& variant : catalog.variants) {
    WriteString(&output, variant.source_id);
    WriteString(&output, variant.allele0);
    WriteString(&output, variant.allele1);
    WriteU64(&output, variant.weights.size());
    for (const auto& weight : variant.weights) {
      WriteU32(&output, weight.score_idx);
      WriteU8(&output, weight.effect_allele_idx);
      WriteDouble(&output, weight.weight);
      ++weights_written;
    }
    ++variants_written;
    if (progress && !(variants_written % 1000000)) {
      const std::streampos position = output.tellp();
      progress->MaybeEvent(
          "compile", "serialize",
          {{"scores_written", catalog.scores.size()},
           {"variants_written", variants_written},
           {"variants_total", catalog.variants.size()},
           {"weights_written", weights_written},
           {"weights_total", catalog.weight_ct},
           {"output_bytes", position >= 0 ? static_cast<uint64_t>(position) : 0}});
    }
  }
  output.close();
  if (!output) {
    throw std::runtime_error("cannot finish " + temporary);
  }
  std::filesystem::rename(temporary, path);
  if (progress) {
    progress->Event(
        "compile", "serialization_complete",
        {{"scores_written", catalog.scores.size()},
         {"variants_written", variants_written},
         {"weights_written", weights_written},
         {"output_bytes", std::filesystem::file_size(path)}});
  }
}

CompiledCatalog ReadCompiledCatalog(const std::string& path,
                                    ProgressReporter* progress) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open " + path);
  }
  char magic[sizeof(kMagic)];
  ReadBytes(&input, magic, sizeof(magic), path);
  if (std::memcmp(magic, kMagic, sizeof(kMagic)) != 0) {
    throw std::runtime_error(path + " is not a pgensparsescore catalog");
  }
  const uint32_t format_version = ReadU32(&input, path);
  if (format_version != 1 && format_version != kFormatVersion) {
    throw std::runtime_error(path + " has an unsupported catalog version");
  }
  const uint32_t score_ct = ReadU32(&input, path);
  const uint64_t variant_ct = ReadU64(&input, path);
  const uint64_t expected_weight_ct = ReadU64(&input, path);
  if (!score_ct || score_ct > std::numeric_limits<uint32_t>::max() ||
      variant_ct > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error(path + " has invalid catalog dimensions");
  }

  CompiledCatalog result;
  result.scores.reserve(score_ct);
  result.variants.reserve(static_cast<size_t>(variant_ct));
  for (uint32_t idx = 0; idx < score_ct; ++idx) {
    ScoreInfo score;
    score.id = ReadString(&input, path);
    score.path = ReadString(&input, path);
    score.input_weight_ct = ReadU64(&input, path);
    score.zero_weight_ct = ReadU64(&input, path);
    score.excluded_weight_ct = ReadU64(&input, path);
    if (format_version >= 2) {
      score.duplicate_weight_ct = ReadU64(&input, path);
    }
    score.catalog_weight_ct = ReadU64(&input, path);
    result.scores.push_back(std::move(score));
  }
  for (uint64_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
    CompiledVariant variant;
    variant.source_id = ReadString(&input, path);
    variant.allele0 = ReadString(&input, path);
    variant.allele1 = ReadString(&input, path);
    const uint64_t weight_ct = ReadU64(&input, path);
    if (weight_ct > expected_weight_ct) {
      throw std::runtime_error(path + " has an invalid variant weight count");
    }
    variant.weights.reserve(static_cast<size_t>(weight_ct));
    uint32_t previous_score_idx = 0;
    for (uint64_t weight_idx = 0; weight_idx < weight_ct; ++weight_idx) {
      const uint32_t score_idx = ReadU32(&input, path);
      const uint8_t effect_allele_idx = ReadU8(&input, path);
      const double weight = ReadDouble(&input, path);
      if (score_idx >= score_ct || effect_allele_idx > 1 ||
          (weight_idx && score_idx < previous_score_idx)) {
        throw std::runtime_error(path + " contains an invalid weight record");
      }
      previous_score_idx = score_idx;
      variant.weights.push_back({score_idx, weight, effect_allele_idx});
      ++result.weight_ct;
    }
    result.variants.push_back(std::move(variant));
    if (progress && !((variant_idx + 1) % 1000000)) {
      progress->MaybeEvent(
          "score", "read_catalog",
          {{"scores_read", score_ct},
           {"variants_read", variant_idx + 1},
           {"variants_total", variant_ct},
           {"weights_read", result.weight_ct},
           {"weights_total", expected_weight_ct}});
    }
  }
  if (result.weight_ct != expected_weight_ct) {
    throw std::runtime_error(path + " weight count disagrees with its header");
  }
  char extra = 0;
  if (input.read(&extra, 1)) {
    throw std::runtime_error(path + " has trailing data");
  }
  if (progress) {
    progress->Event("score", "catalog_loaded",
                    {{"scores_read", score_ct},
                     {"variants_read", variant_ct},
                     {"weights_read", result.weight_ct}});
  }
  return result;
}

Catalog MaterializeCompiledCatalog(const CompiledCatalog& compiled,
                                   const std::vector<Variant>& variants,
                                   const VariantMap* variant_map,
                                   ProgressReporter* progress) {
  std::unordered_map<std::string, uint32_t> variant_by_id;
  variant_by_id.reserve(variants.size());
  for (uint32_t idx = 0; idx < variants.size(); ++idx) {
    if (variants[idx].id.empty() || variants[idx].id == ".") {
      continue;
    }
    if (!variant_by_id.emplace(variants[idx].id, idx).second) {
      throw std::runtime_error("duplicate nonmissing PVAR ID: " +
                               variants[idx].id);
    }
  }

  Catalog result;
  result.scores = compiled.scores;
  result.intercepts.assign(result.scores.size(), 0.0);
  for (auto& score : result.scores) {
    score.matched_weight_ct = 0;
    score.missing_variant_ct = 0;
    score.missing_frequency_ct = 0;
    score.ref_effect_ct = 0;
    score.alt_effect_ct = 0;
    score.ref_effect_intercept = 0.0;
  }

  uint64_t source_variants_processed = 0;
  uint64_t matched_weights = 0;
  for (const auto& source_variant : compiled.variants) {
    ++source_variants_processed;
    if (progress && !(source_variants_processed % 1000000)) {
      progress->MaybeEvent(
          "score", "materialize_catalog",
          {{"source_variants_processed", source_variants_processed},
           {"source_variants_total", compiled.variants.size()},
           {"matched_variants", result.variants.size()},
           {"matched_weights", matched_weights}});
    }
    std::string target_id = source_variant.source_id;
    if (variant_map) {
      const auto mapping = variant_map->find(source_variant.source_id);
      if (mapping == variant_map->end()) {
        for (const auto& weight : source_variant.weights) {
          ++result.scores[weight.score_idx].missing_variant_ct;
        }
        continue;
      }
      target_id = mapping->second;
    }
    const auto pvar = variant_by_id.find(target_id);
    if (pvar == variant_by_id.end()) {
      for (const auto& weight : source_variant.weights) {
        ++result.scores[weight.score_idx].missing_variant_ct;
      }
      continue;
    }
    const Variant& target = variants[pvar->second];
    if (target.alt.find(',') != std::string::npos) {
      throw std::runtime_error("compiled catalog references multiallelic PVAR " +
                               target.id);
    }
    const bool allele_pair_matches =
        (target.ref == source_variant.allele0 &&
         target.alt == source_variant.allele1) ||
        (target.ref == source_variant.allele1 &&
         target.alt == source_variant.allele0);
    if (!allele_pair_matches) {
      throw std::runtime_error(
          "compiled catalog alleles disagree with PVAR for " + target.id +
          " (PVAR " + target.ref + "/" + target.alt + ", catalog " +
          source_variant.allele0 + "/" + source_variant.allele1 + ")");
    }

    VariantEdges materialized;
    materialized.variant_idx = pvar->second;
    materialized.edges.reserve(source_variant.weights.size());
    for (const auto& weight : source_variant.weights) {
      const std::string& effect = weight.effect_allele_idx
                                      ? source_variant.allele1
                                      : source_variant.allele0;
      Edge edge;
      edge.score_idx = weight.score_idx;
      if (effect == target.alt) {
        edge.beta_alt = weight.weight;
        ++result.scores[weight.score_idx].alt_effect_ct;
      } else {
        edge.beta_alt = -weight.weight;
        edge.ref_effect = true;
        result.intercepts[weight.score_idx] += 2.0 * weight.weight;
        ++result.scores[weight.score_idx].ref_effect_ct;
        result.scores[weight.score_idx].ref_effect_intercept +=
            2.0 * weight.weight;
      }
      ++result.scores[weight.score_idx].matched_weight_ct;
      ++matched_weights;
      materialized.edges.push_back(edge);
    }
    result.variants.push_back(std::move(materialized));
  }
  std::sort(result.variants.begin(), result.variants.end(),
            [](const VariantEdges& lhs, const VariantEdges& rhs) {
              return lhs.variant_idx < rhs.variant_idx;
            });
  if (progress) {
    progress->Event("score", "catalog_materialized",
                    {{"source_variants_processed", compiled.variants.size()},
                     {"matched_variants", result.variants.size()},
                     {"matched_weights", matched_weights}});
  }
  return result;
}

}  // namespace pgensparsescore
