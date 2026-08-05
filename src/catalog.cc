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
#include "compiled_catalog.h"

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

}  // namespace

PvarData ReadPvar(
    const std::string& path,
    const std::unordered_set<std::string>* included_ids) {
  LineReader reader(path);
  std::string line;
  Header header;
  PvarData result;
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
    if (result.row_ct == std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error(path + " exceeds the supported PVAR row count");
    }
    const uint32_t pgen_variant_idx = result.row_ct++;
    if (included_ids && !included_ids->count(fields[id_idx])) {
      continue;
    }
    result.variants.push_back({fields[chrom_idx], fields[id_idx],
                               fields[ref_idx], fields[alt_idx],
                               pgen_variant_idx});
  }
  if (header.empty()) {
    throw std::runtime_error(path + " has no header");
  }
  if (!result.row_ct) {
    throw std::runtime_error(path + " has no variants");
  }
  return result;
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
    samples.push_back(
        {fid_iter == header.end()
             ? std::nullopt
             : std::optional<std::string>(fields[fid_iter->second]),
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

VariantMap ReadVariantMap(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) {
    throw std::runtime_error(path + " is empty");
  }
  const Header header = MakeHeader(SplitTabs(line));
  const size_t source_idx = RequireColumn(header, "SOURCE_ID", path);
  const size_t target_idx = RequireColumn(header, "TARGET_ID", path);
  VariantMap result;
  std::unordered_set<std::string> targets;
  uint64_t line_number = 1;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty()) {
      continue;
    }
    const auto fields = SplitTabs(line);
    RequireFieldCount(fields, std::max(source_idx, target_idx), path,
                      line_number);
    const std::string& source = fields[source_idx];
    const std::string& target = fields[target_idx];
    if (source.empty() || target.empty()) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has a blank variant ID");
    }
    if (!result.emplace(source, target).second) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " duplicates SOURCE_ID " + source);
    }
    if (!targets.insert(target).second) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " duplicates TARGET_ID " + target);
    }
  }
  if (result.empty()) {
    throw std::runtime_error(path + " has no variant mappings");
  }
  return result;
}

Catalog CompileCatalog(const std::string& manifest_path,
                       const std::vector<Variant>& variants,
                       const VariantMap* variant_map) {
  const auto compiled = CompileSourceCatalog(manifest_path);
  return MaterializeCompiledCatalog(compiled, variants, variant_map);
}

}  // namespace pgensparsescore
