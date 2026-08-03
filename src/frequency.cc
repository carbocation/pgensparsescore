// SPDX-License-Identifier: GPL-3.0-only
#include "frequency.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "io.h"

namespace pgensparsescore {

namespace {

using Header = std::unordered_map<std::string, size_t>;

std::string Upper(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return value;
}

Header MakeHeader(const std::vector<std::string>& fields,
                  const std::string& path) {
  Header result;
  for (size_t idx = 0; idx < fields.size(); ++idx) {
    std::string name = fields[idx];
    if (!name.empty() && name.front() == '#') {
      name.erase(0, 1);
    }
    name = Upper(name);
    if (!result.emplace(name, idx).second) {
      throw std::runtime_error(path + " has duplicate frequency column " + name);
    }
  }
  return result;
}

size_t Require(const Header& header, const std::string& name,
               const std::string& path) {
  const auto iter = header.find(name);
  if (iter == header.end()) {
    throw std::runtime_error(path + " is missing frequency column " + name);
  }
  return iter->second;
}

size_t Find(const Header& header,
            std::initializer_list<const char*> candidates) {
  for (const char* candidate : candidates) {
    const auto iter = header.find(candidate);
    if (iter != header.end()) {
      return iter->second;
    }
  }
  return static_cast<size_t>(-1);
}

double ParseNumber(const std::string& value, const std::string& path,
                   uint64_t line_number, const std::string& column) {
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(value.c_str(), &end);
  if (errno || end == value.c_str() || *end != '\0' || !std::isfinite(parsed)) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has invalid " + column + ": " + value);
  }
  return parsed;
}

double ParseBiallelicValue(const std::string& raw, const std::string& alt,
                           const std::string& path, uint64_t line_number,
                           const std::string& column) {
  if (raw.find(',') != std::string::npos) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has multiallelic " + column + ": " + raw);
  }
  const size_t equals = raw.find('=');
  if (equals == std::string::npos) {
    return ParseNumber(raw, path, line_number, column);
  }
  if (Upper(raw.substr(0, equals)) != alt) {
    throw std::runtime_error(path + ": line " + std::to_string(line_number) +
                             " has " + column + " for the wrong ALT allele: " +
                             raw);
  }
  return ParseNumber(raw.substr(equals + 1), path, line_number, column);
}

}  // namespace

FrequencyTable ReadFrequencyTable(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) {
    throw std::runtime_error(path + " is empty");
  }
  const Header header = MakeHeader(SplitTabs(line), path);
  const size_t id_idx = Require(header, "ID", path);
  const size_t ref_idx = Require(header, "REF", path);
  const size_t alt_idx = Find(header, {"ALT", "ALT1"});
  if (alt_idx == static_cast<size_t>(-1)) {
    throw std::runtime_error(path + " is missing frequency column ALT or ALT1");
  }
  const size_t mean_idx = Find(header, {"ALT_DOSAGE_MEAN"});
  const size_t freq_idx = Find(header, {"ALT_FREQS", "ALT1_FREQ"});
  const size_t alt_count_idx = Find(header, {"ALT_CTS", "ALT1_CT"});
  const size_t ref_count_idx = Find(header, {"REF_CT"});
  const size_t obs_idx = Find(header, {"OBS_CT"});
  if (mean_idx == static_cast<size_t>(-1) &&
      freq_idx == static_cast<size_t>(-1) &&
      alt_count_idx == static_cast<size_t>(-1) &&
      ref_count_idx == static_cast<size_t>(-1)) {
    throw std::runtime_error(
        path + " has no ALT_DOSAGE_MEAN, ALT_FREQS, ALT_CTS, or REF_CT column");
  }
  if (mean_idx == static_cast<size_t>(-1) &&
      freq_idx == static_cast<size_t>(-1) &&
      obs_idx == static_cast<size_t>(-1)) {
    throw std::runtime_error(path +
                             " requires OBS_CT when allele counts are used");
  }

  const size_t max_idx = std::max(
      {id_idx, ref_idx, alt_idx,
       mean_idx == static_cast<size_t>(-1) ? 0 : mean_idx,
       freq_idx == static_cast<size_t>(-1) ? 0 : freq_idx,
       alt_count_idx == static_cast<size_t>(-1) ? 0 : alt_count_idx,
       ref_count_idx == static_cast<size_t>(-1) ? 0 : ref_count_idx,
       obs_idx == static_cast<size_t>(-1) ? 0 : obs_idx});

  FrequencyTable result;
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
    const std::string ref = Upper(fields[ref_idx]);
    const std::string alt = Upper(fields[alt_idx]);
    if (ref.empty() || alt.empty() || alt.find(',') != std::string::npos) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " is not biallelic");
    }

    double alt_dosage_mean = 0.0;
    if (mean_idx != static_cast<size_t>(-1)) {
      alt_dosage_mean = ParseNumber(fields[mean_idx], path, line_number,
                                    "ALT_DOSAGE_MEAN");
    } else if (freq_idx != static_cast<size_t>(-1)) {
      alt_dosage_mean = 2.0 * ParseBiallelicValue(
          fields[freq_idx], alt, path, line_number, "ALT_FREQS");
    } else {
      const double observations =
          ParseNumber(fields[obs_idx], path, line_number, "OBS_CT");
      if (observations <= 0.0) {
        throw std::runtime_error(path + ": line " +
                                 std::to_string(line_number) +
                                 " has nonpositive OBS_CT");
      }
      const double alt_count =
          alt_count_idx != static_cast<size_t>(-1)
              ? ParseBiallelicValue(fields[alt_count_idx], alt, path,
                                    line_number, "ALT_CTS")
              : observations - ParseNumber(fields[ref_count_idx], path,
                                            line_number, "REF_CT");
      alt_dosage_mean = 2.0 * alt_count / observations;
    }
    if (alt_dosage_mean < 0.0 || alt_dosage_mean > 2.0) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has ALT dosage mean outside [0,2]");
    }
    if (!result.emplace(fields[id_idx],
                        AlleleFrequency{ref, alt, alt_dosage_mean})
             .second) {
      throw std::runtime_error(path + ": duplicate frequency ID " +
                               fields[id_idx]);
    }
  }
  if (result.empty()) {
    throw std::runtime_error(path + " has no frequency rows");
  }
  return result;
}

}  // namespace pgensparsescore
