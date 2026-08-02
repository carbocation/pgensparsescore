// SPDX-License-Identifier: GPL-3.0-only
#include "catalog.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "io.h"

namespace pgensparsescore {

namespace {

using Header = std::unordered_map<std::string, size_t>;

Header MakeHeader(const std::vector<std::string>& fields) {
  Header header;
  for (size_t idx = 0; idx < fields.size(); ++idx) {
    std::string name = fields[idx];
    if (!name.empty() && name[0] == '#') {
      name.erase(0, 1);
    }
    if (!header.emplace(name, idx).second) {
      throw std::runtime_error("duplicate column in header: " + name);
    }
  }
  return header;
}

size_t RequireColumn(const Header& header, const std::string& name,
                     const std::string& path) {
  const auto iter = header.find(name);
  if (iter == header.end()) {
    throw std::runtime_error(path + " is missing required column " + name);
  }
  return iter->second;
}

void RequireFieldCount(const std::vector<std::string>& fields,
                       size_t required_idx, const std::string& path,
                       uint64_t line_number) {
  if (fields.size() <= required_idx) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has too few fields");
  }
}

double ParseDouble(const std::string& value, const std::string& path,
                   uint64_t line_number) {
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(value.c_str(), &end);
  if (errno || end == value.c_str() || *end != '\0' || !std::isfinite(parsed)) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has invalid finite weight: " + value);
  }
  return parsed;
}

std::string FindFirst(const Header& header,
                      std::initializer_list<const char*> candidates,
                      const std::string& path) {
  for (const char* candidate : candidates) {
    if (header.count(candidate)) {
      return candidate;
    }
  }
  std::string message = path + " is missing required column (one of";
  for (const char* candidate : candidates) {
    message += " ";
    message += candidate;
  }
  message += ")";
  throw std::runtime_error(message);
}

}  // namespace

std::vector<Variant> ReadPvar(const std::string& path) {
  LineReader reader(path);
  std::string line;
  Header header;
  std::vector<Variant> variants;
  uint64_t line_number = 0;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty() || line.rfind("##", 0) == 0) {
      continue;
    }
    if (header.empty()) {
      header = MakeHeader(SplitTabs(line));
      continue;
    }
    const auto fields = SplitTabs(line);
    const size_t chrom_idx = RequireColumn(header, "CHROM", path);
    const size_t id_idx = RequireColumn(header, "ID", path);
    const size_t ref_idx = RequireColumn(header, "REF", path);
    const size_t alt_idx = RequireColumn(header, "ALT", path);
    const size_t max_idx = std::max({chrom_idx, id_idx, ref_idx, alt_idx});
    RequireFieldCount(fields, max_idx, path, line_number);
    variants.push_back({fields[chrom_idx], fields[id_idx], fields[ref_idx],
                        fields[alt_idx]});
  }
  if (header.empty()) {
    throw std::runtime_error(path + " has no header");
  }
  if (variants.empty()) {
    throw std::runtime_error(path + " has no variants");
  }
  return variants;
}

std::vector<Sample> ReadPsam(const std::string& path) {
  LineReader reader(path);
  std::string line;
  Header header;
  std::vector<Sample> samples;
  uint64_t line_number = 0;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty() || line.rfind("##", 0) == 0) {
      continue;
    }
    if (header.empty()) {
      header = MakeHeader(SplitTabs(line));
      continue;
    }
    const auto fields = SplitTabs(line);
    const size_t iid_idx = RequireColumn(header, "IID", path);
    const auto fid_iter = header.find("FID");
    const size_t max_idx = fid_iter == header.end()
                               ? iid_idx
                               : std::max(iid_idx, fid_iter->second);
    RequireFieldCount(fields, max_idx, path, line_number);
    samples.push_back({fid_iter == header.end() ? "0" : fields[fid_iter->second],
                       fields[iid_idx]});
  }
  if (header.empty()) {
    throw std::runtime_error(path + " has no header");
  }
  if (samples.empty()) {
    throw std::runtime_error(path + " has no samples");
  }
  return samples;
}

Catalog CompileCatalog(const std::string& manifest_path,
                       const std::vector<Variant>& variants) {
  std::unordered_map<std::string, uint32_t> variant_by_id;
  variant_by_id.reserve(variants.size());
  for (uint32_t idx = 0; idx < variants.size(); ++idx) {
    const auto& variant = variants[idx];
    if (variant.id.empty() || variant.id == ".") {
      continue;
    }
    if (!variant_by_id.emplace(variant.id, idx).second) {
      throw std::runtime_error("duplicate nonmissing PVAR ID: " + variant.id);
    }
  }

  struct ManifestRow {
    std::string score;
    std::string path;
  };
  std::vector<ManifestRow> manifest;
  {
    LineReader reader(manifest_path);
    std::string line;
    if (!reader.GetLine(&line)) {
      throw std::runtime_error(manifest_path + " is empty");
    }
    const Header header = MakeHeader(SplitTabs(line));
    const std::string score_col =
        FindFirst(header, {"SCORE", "SCORE_ID"}, manifest_path);
    const size_t score_idx = header.at(score_col);
    const size_t path_idx = RequireColumn(header, "PATH", manifest_path);
    uint64_t line_number = 1;
    std::unordered_set<std::string> seen_scores;
    const std::filesystem::path base =
        std::filesystem::absolute(manifest_path).parent_path();
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) {
        continue;
      }
      const auto fields = SplitTabs(line);
      RequireFieldCount(fields, std::max(score_idx, path_idx), manifest_path,
                        line_number);
      if (!seen_scores.insert(fields[score_idx]).second) {
        throw std::runtime_error("duplicate manifest SCORE: " +
                                 fields[score_idx]);
      }
      std::filesystem::path weight_path(fields[path_idx]);
      if (weight_path.is_relative()) {
        weight_path = base / weight_path;
      }
      manifest.push_back({fields[score_idx], weight_path.lexically_normal().string()});
    }
  }
  if (manifest.empty()) {
    throw std::runtime_error(manifest_path + " has no scores");
  }

  Catalog catalog;
  catalog.scores.reserve(manifest.size());
  catalog.intercepts.assign(manifest.size(), 0.0);
  std::unordered_map<uint32_t, std::vector<Edge>> edges_by_variant;

  for (uint32_t score_idx = 0; score_idx < manifest.size(); ++score_idx) {
    const auto& item = manifest[score_idx];
    ScoreInfo info;
    info.id = item.score;
    info.path = item.path;
    LineReader reader(item.path);
    std::string line;
    if (!reader.GetLine(&line)) {
      throw std::runtime_error(item.path + " is empty");
    }
    const Header header = MakeHeader(SplitTabs(line));
    const size_t snp_idx = RequireColumn(header, "SNP", item.path);
    const size_t effect_idx = RequireColumn(header, "EFFECT_ALLELE", item.path);
    const size_t other_idx = RequireColumn(header, "OTHER_ALLELE", item.path);
    const size_t weight_idx =
        RequireColumn(header, "EFFECT_ALLELE_WEIGHT", item.path);
    const size_t max_idx = std::max({snp_idx, effect_idx, other_idx, weight_idx});
    std::unordered_set<uint32_t> score_variants;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) {
        continue;
      }
      ++info.input_weight_ct;
      const auto fields = SplitTabs(line);
      RequireFieldCount(fields, max_idx, item.path, line_number);
      const auto pvar_iter = variant_by_id.find(fields[snp_idx]);
      if (pvar_iter == variant_by_id.end()) {
        ++info.missing_variant_ct;
        continue;
      }
      const uint32_t variant_idx = pvar_iter->second;
      if (!score_variants.insert(variant_idx).second) {
        throw std::runtime_error(item.path + ": line " +
                                 std::to_string(line_number) +
                                 " duplicates PVAR variant " + fields[snp_idx]);
      }
      const Variant& variant = variants[variant_idx];
      if (variant.alt.find(',') != std::string::npos) {
        throw std::runtime_error(item.path + ": line " +
                                 std::to_string(line_number) +
                                 " references multiallelic PVAR variant " +
                                 variant.id);
      }
      const std::string& effect = fields[effect_idx];
      const std::string& other = fields[other_idx];
      const double weight = ParseDouble(fields[weight_idx], item.path, line_number);
      double beta_alt = 0.0;
      if (effect == variant.alt && other == variant.ref) {
        beta_alt = weight;
        ++info.alt_effect_ct;
      } else if (effect == variant.ref && other == variant.alt) {
        beta_alt = -weight;
        info.ref_effect_intercept += 2.0 * weight;
        ++info.ref_effect_ct;
      } else {
        throw std::runtime_error(
            item.path + ": line " + std::to_string(line_number) +
            " allele mismatch for " + variant.id + " (PVAR " + variant.ref +
            "/" + variant.alt + ", weight " + effect + "/" + other + ")");
      }
      edges_by_variant[variant_idx].push_back({score_idx, beta_alt});
      ++info.matched_weight_ct;
    }
    catalog.intercepts[score_idx] = info.ref_effect_intercept;
    catalog.scores.push_back(std::move(info));
  }

  catalog.variants.reserve(edges_by_variant.size());
  for (auto& [variant_idx, edges] : edges_by_variant) {
    std::sort(edges.begin(), edges.end(),
              [](const Edge& lhs, const Edge& rhs) {
                return lhs.score_idx < rhs.score_idx;
              });
    catalog.variants.push_back({variant_idx, std::move(edges)});
  }
  std::sort(catalog.variants.begin(), catalog.variants.end(),
            [](const VariantEdges& lhs, const VariantEdges& rhs) {
              return lhs.variant_idx < rhs.variant_idx;
            });
  return catalog;
}

}  // namespace pgensparsescore
