// SPDX-License-Identifier: GPL-3.0-only
#include "pfile.h"

#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <unordered_map>

#include "io.h"

namespace pgensparsescore {

std::vector<PfileSpec> ReadPfileList(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) {
    throw std::runtime_error(path + " is empty");
  }
  const auto header_fields = SplitTabs(line);
  std::unordered_map<std::string, size_t> header;
  for (size_t idx = 0; idx < header_fields.size(); ++idx) {
    std::string name = header_fields[idx];
    if (!name.empty() && name.front() == '#') {
      name.erase(0, 1);
    }
    if (!header.emplace(name, idx).second) {
      throw std::runtime_error(path + " has duplicate column " + name);
    }
  }
  for (const char* name : {"PGEN", "PVAR", "PSAM"}) {
    if (!header.count(name)) {
      throw std::runtime_error(path + " is missing column " + name);
    }
  }
  const size_t max_idx =
      std::max({header.at("PGEN"), header.at("PVAR"), header.at("PSAM")});
  const std::filesystem::path base =
      std::filesystem::absolute(path).parent_path();
  auto resolve = [&](const std::string& value) {
    std::filesystem::path resolved(value);
    if (resolved.is_relative()) {
      resolved = base / resolved;
    }
    return resolved.lexically_normal().string();
  };

  std::vector<PfileSpec> result;
  uint64_t line_number = 1;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty()) {
      continue;
    }
    const auto fields = SplitTabs(line);
    if (fields.size() <= max_idx) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has too few fields");
    }
    result.push_back({resolve(fields[header.at("PGEN")]),
                      resolve(fields[header.at("PVAR")]),
                      resolve(fields[header.at("PSAM")])});
  }
  if (result.empty()) {
    throw std::runtime_error(path + " has no PGEN rows");
  }
  return result;
}

}  // namespace pgensparsescore
