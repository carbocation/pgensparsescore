// SPDX-License-Identifier: GPL-3.0-only
#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <sched.h>
#endif

#include "catalog.h"
#include "compiled_catalog.h"
#include "frequency.h"
#include "fragment_scorer.h"
#include "io.h"
#include "mapped_matrix.h"
#include "pfile.h"
#include "pgen_reader.h"
#include "progress.h"
#include "scorer.h"
#include "score_fragment.h"
#include "support_index.h"
#include "variant_index.h"

namespace {

enum class ScoreOutputFormat {
  kWideTsv,
  kScoreMajorBinary,
};

struct Options {
  std::string pgen;
  std::string pvar;
  std::string psam;
  std::string pfile_list;
  std::string manifest;
  std::string compiled_catalog;
  std::string fragment_list;
  std::string variant_index;
  std::string support_index;
  std::string score_schema;
  std::string variant_map;
  std::string read_freq;
  std::string progress_jsonl;
  std::string progress_interval_seconds;
  std::string threads;
  std::string dense_kernel;
  std::string output_format;
  std::string out;
  pgensparsescore::MissingFrequencyPolicy missing_frequency_policy =
      pgensparsescore::MissingFrequencyPolicy::kCohort;
  bool missing_frequency_policy_supplied = false;
};

struct CompileOptions {
  std::string manifest;
  std::string variant_map;
  std::string progress_jsonl;
  std::string progress_interval_seconds;
  std::string out;
};

struct VariantIndexOptions {
  pgensparsescore::VariantIndexBuildOptions build;
  std::string progress_jsonl;
  std::string progress_interval_seconds;
};

struct FragmentCompileOptions {
  pgensparsescore::ScoreFragmentCompileOptions build;
  std::string progress_jsonl;
  std::string progress_interval_seconds;
};

struct VariantBitsMergeOptions {
  std::string list_path;
  std::string output_path;
};

struct SupportIndexOptions {
  pgensparsescore::SupportIndexBuildOptions build;
  std::string progress_jsonl;
  std::string progress_interval_seconds;
};

struct FragmentSupportOptions {
  std::string fragment_path;
  std::string support_index_path;
  std::string output_path;
};

class RemoveFileOnExit {
 public:
  explicit RemoveFileOnExit(std::string path) : path_(std::move(path)) {}
  ~RemoveFileOnExit() {
    if (active_) {
      std::error_code ignored;
      std::filesystem::remove(path_, ignored);
    }
  }

  void RemoveNow() {
    if (!std::filesystem::remove(path_)) {
      throw std::runtime_error("cannot remove temporary file " + path_);
    }
    active_ = false;
  }

  void Release() { active_ = false; }

 private:
  std::string path_;
  bool active_ = true;
};

void PrintUsage(std::ostream& stream) {
  stream <<
      "Usage:\n"
      "  pgensparsescore build-variant-index --variant-list FILE \\\n"
      "                          [--source-id-column NAME] \\\n"
      "                          [--target-id-column NAME] \\\n"
      "                          [--strip-target-id-prefix TEXT] \\\n"
      "                          [--ref-column NAME] [--alt-column NAME] \\\n"
      "                          [--block-size N] [--progress-jsonl FILE] \\\n"
      "                          [--progress-interval-seconds N] \\\n"
      "                          --out VARIANTS.index.bin\n"
      "  pgensparsescore build-support-index --variant-index FILE \\\n"
      "                          (--pvar FILE | --pvar-list FILE) \\\n"
      "                          --read-freq FILE \\\n"
      "                          [--progress-jsonl FILE] \\\n"
      "                          [--progress-interval-seconds N] \\\n"
      "                          --out VARIANTS.support.bin\n"
      "  pgensparsescore compile-fragment --manifest FILE \\\n"
      "                          --variant-index FILE \\\n"
      "                          [--support-index FILE] [--temp-dir DIR] \\\n"
      "                          [--minimum-supported-fraction X] \\\n"
      "                          [--progress-jsonl FILE] \\\n"
      "                          [--progress-interval-seconds N] \\\n"
      "                          --out SCORES.fragment.bin\n"
      "  pgensparsescore merge-variant-bits --list FILE \\\n"
      "                          --out VARIANTS.bits\n"
      "  pgensparsescore report-fragment-support --fragment FILE \\\n"
      "                          --support-index FILE \\\n"
      "                          --out REPORT.tsv\n"
      "  pgensparsescore compile --manifest FILE [--variant-map FILE] \\\n"
      "                          [--progress-jsonl FILE] \\\n"
      "                          [--progress-interval-seconds N] \\\n"
      "                          --out CATALOG.bin\n"
      "  pgensparsescore --pgen FILE --pvar FILE --psam FILE \\\n"
      "                     (--manifest FILE | --compiled-catalog FILE) \\\n"
      "                     [--variant-map FILE] \\\n"
      "                     [--read-freq FILE] \\\n"
      "                     [--progress-jsonl FILE] \\\n"
      "                     [--progress-interval-seconds N] \\\n"
      "                     [--output-format wide-tsv|score-major-bin] \\\n"
      "                     [--missing-freq cohort|error|omit] --out PREFIX\n"
      "  pgensparsescore --pfile-list FILE \\\n"
      "                     (--manifest FILE | --compiled-catalog FILE) \\\n"
      "                     [--variant-map FILE] \\\n"
      "                     [--read-freq FILE] \\\n"
      "                     [--progress-jsonl FILE] \\\n"
      "                     [--progress-interval-seconds N] \\\n"
      "                     [--output-format wide-tsv|score-major-bin] \\\n"
      "                     [--missing-freq cohort|error|omit] --out PREFIX\n";
  stream <<
      "  pgensparsescore (--pgen FILE --pvar FILE --psam FILE | \\\n"
      "                     --pfile-list FILE) --fragment-list FILE \\\n"
      "                     --variant-index FILE [--read-freq FILE] \\\n"
      "                     [--support-index FILE] \\\n"
      "                     [--score-schema FILE] \\\n"
      "                     [--progress-jsonl FILE] \\\n"
      "                     [--progress-interval-seconds N] \\\n"
      "                     [--threads N] \\\n"
      "                     [--dense-kernel auto|direct|onemkl] \\\n"
      "                     [--output-format wide-tsv|score-major-bin] \\\n"
      "                     [--missing-freq cohort|error|omit] --out PREFIX\n";
}

ScoreOutputFormat ParseScoreOutputFormat(const std::string& value) {
  if (value.empty() || value == "wide-tsv") {
    return ScoreOutputFormat::kWideTsv;
  }
  if (value == "score-major-bin") {
    return ScoreOutputFormat::kScoreMajorBinary;
  }
  throw std::runtime_error(
      "--output-format must be wide-tsv or score-major-bin");
}

const char* ScoreOutputFormatName(ScoreOutputFormat format) {
  switch (format) {
    case ScoreOutputFormat::kWideTsv:
      return "pgensparsescore-wide-tsv-v1";
    case ScoreOutputFormat::kScoreMajorBinary:
      return "pgensparsescore-score-major-f64-v1";
  }
  throw std::runtime_error("invalid score output format");
}

std::string ScoreOutputPath(const std::string& prefix,
                            ScoreOutputFormat format) {
  return prefix + (format == ScoreOutputFormat::kWideTsv
                       ? ".scores.tsv.gz"
                       : ".scores.f64le");
}

pgensparsescore::DenseScoringKernel ParseDenseScoringKernel(
    const std::string& value) {
  if (value.empty() || value == "auto") {
    return pgensparsescore::DenseScoringKernel::kAuto;
  }
  if (value == "direct") {
    return pgensparsescore::DenseScoringKernel::kDirect;
  }
  if (value == "onemkl") {
    return pgensparsescore::DenseScoringKernel::kOneMkl;
  }
  throw std::runtime_error(
      "--dense-kernel must be auto, direct, or onemkl");
}

uint32_t ParsePositiveU32(const std::string& value,
                          const std::string& argument) {
  uint32_t result = 0;
  const auto parsed =
      std::from_chars(value.data(), value.data() + value.size(), result);
  if (parsed.ec != std::errc() || parsed.ptr != value.data() + value.size() ||
      !result) {
    throw std::runtime_error(argument + " must be a positive integer");
  }
  return result;
}

double ParseFraction(const std::string& value, const std::string& argument) {
  char* end = nullptr;
  errno = 0;
  const double result = std::strtod(value.c_str(), &end);
  if (errno || end != value.c_str() + value.size() || !std::isfinite(result) ||
      result < 0.0 || result > 1.0) {
    throw std::runtime_error(argument + " must be a number from 0 through 1");
  }
  return result;
}

uint32_t PhysicalCoreCount() {
#if defined(__APPLE__)
  uint32_t physical_core_ct = 0;
  size_t byte_ct = sizeof(physical_core_ct);
  if (!sysctlbyname("hw.physicalcpu", &physical_core_ct, &byte_ct, nullptr,
                    0) &&
      physical_core_ct) {
    return physical_core_ct;
  }
#elif defined(__linux__)
  cpu_set_t allowed_cpus;
  CPU_ZERO(&allowed_cpus);
  if (!sched_getaffinity(0, sizeof(allowed_cpus), &allowed_cpus)) {
    std::set<std::pair<int, int>> physical_cores;
    for (int cpu_idx = 0; cpu_idx < CPU_SETSIZE; ++cpu_idx) {
      if (!CPU_ISSET(cpu_idx, &allowed_cpus)) continue;
      const std::string topology = "/sys/devices/system/cpu/cpu" +
                                   std::to_string(cpu_idx) + "/topology/";
      std::ifstream package_input(topology + "physical_package_id");
      std::ifstream core_input(topology + "core_id");
      int package_idx = 0;
      int core_idx = 0;
      if (package_input >> package_idx && core_input >> core_idx) {
        physical_cores.emplace(package_idx, core_idx);
      }
    }
    if (!physical_cores.empty()) return physical_cores.size();
  }
#endif
  const uint32_t logical_core_ct = std::thread::hardware_concurrency();
  return logical_core_ct ? logical_core_ct : 1;
}

uint32_t LogicalCoreCount() {
#ifdef __linux__
  cpu_set_t allowed_cpus;
  CPU_ZERO(&allowed_cpus);
  if (!sched_getaffinity(0, sizeof(allowed_cpus), &allowed_cpus)) {
    uint32_t allowed_cpu_ct = 0;
    for (int cpu_idx = 0; cpu_idx < CPU_SETSIZE; ++cpu_idx) {
      allowed_cpu_ct += CPU_ISSET(cpu_idx, &allowed_cpus) ? 1 : 0;
    }
    if (allowed_cpu_ct) return allowed_cpu_ct;
  }
#endif
  const uint32_t logical_core_ct = std::thread::hardware_concurrency();
  return logical_core_ct ? logical_core_ct : PhysicalCoreCount();
}

uint32_t ParseProgressInterval(const std::string& value) {
  if (value.empty()) {
    return 30;
  }
  uint32_t result = 0;
  const auto parsed =
      std::from_chars(value.data(), value.data() + value.size(), result);
  if (parsed.ec != std::errc() || parsed.ptr != value.data() + value.size() ||
      !result) {
    throw std::runtime_error(
        "--progress-interval-seconds must be a positive integer");
  }
  return result;
}

VariantIndexOptions ParseVariantIndexOptions(int argc, char** argv) {
  if (argc == 3 && std::string(argv[2]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  VariantIndexOptions options;
  std::string block_size;
  std::unordered_map<std::string, std::string*> destinations{
      {"--variant-list", &options.build.input_path},
      {"--source-id-column", &options.build.source_id_column},
      {"--target-id-column", &options.build.target_id_column},
      {"--strip-target-id-prefix", &options.build.target_id_prefix_to_strip},
      {"--ref-column", &options.build.ref_column},
      {"--alt-column", &options.build.alt_column},
      {"--block-size", &block_size},
      {"--progress-jsonl", &options.progress_jsonl},
      {"--progress-interval-seconds", &options.progress_interval_seconds},
      {"--out", &options.build.output_path},
  };
  std::unordered_set<std::string> seen_arguments;
  for (int idx = 2; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    const auto iter = destinations.find(key);
    if (iter == destinations.end()) {
      throw std::runtime_error("unknown build-variant-index argument: " + key);
    }
    if (++idx >= argc) throw std::runtime_error("missing value after " + key);
    if (!seen_arguments.insert(key).second) {
      throw std::runtime_error("argument supplied twice: " + key);
    }
    *iter->second = argv[idx];
  }
  if (options.build.input_path.empty() || options.build.output_path.empty()) {
    throw std::runtime_error(
        "build-variant-index requires --variant-list and --out");
  }
  if (!block_size.empty()) {
    options.build.block_size = ParsePositiveU32(block_size, "--block-size");
  }
  if (!options.progress_interval_seconds.empty() &&
      options.progress_jsonl.empty()) {
    throw std::runtime_error(
        "--progress-interval-seconds requires --progress-jsonl");
  }
  if (!options.progress_interval_seconds.empty()) {
    ParseProgressInterval(options.progress_interval_seconds);
  }
  return options;
}

FragmentCompileOptions ParseFragmentCompileOptions(int argc, char** argv) {
  if (argc == 3 && std::string(argv[2]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  FragmentCompileOptions options;
  std::string minimum_supported_fraction;
  std::unordered_map<std::string, std::string*> destinations{
      {"--manifest", &options.build.manifest_path},
      {"--variant-index", &options.build.variant_index_path},
      {"--support-index", &options.build.support_index_path},
      {"--minimum-supported-fraction", &minimum_supported_fraction},
      {"--temp-dir", &options.build.temporary_directory},
      {"--progress-jsonl", &options.progress_jsonl},
      {"--progress-interval-seconds", &options.progress_interval_seconds},
      {"--out", &options.build.output_path},
  };
  for (int idx = 2; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    const auto iter = destinations.find(key);
    if (iter == destinations.end()) {
      throw std::runtime_error("unknown compile-fragment argument: " + key);
    }
    if (++idx >= argc) throw std::runtime_error("missing value after " + key);
    if (!iter->second->empty()) {
      throw std::runtime_error("argument supplied twice: " + key);
    }
    *iter->second = argv[idx];
  }
  if (options.build.manifest_path.empty() ||
      options.build.variant_index_path.empty() ||
      options.build.output_path.empty()) {
    throw std::runtime_error(
        "compile-fragment requires --manifest, --variant-index, and --out");
  }
  if (!options.progress_interval_seconds.empty() &&
      options.progress_jsonl.empty()) {
    throw std::runtime_error(
        "--progress-interval-seconds requires --progress-jsonl");
  }
  if (!options.progress_interval_seconds.empty()) {
    ParseProgressInterval(options.progress_interval_seconds);
  }
  if (!minimum_supported_fraction.empty()) {
    options.build.minimum_supported_fraction = ParseFraction(
        minimum_supported_fraction, "--minimum-supported-fraction");
  }
  if (options.build.minimum_supported_fraction > 0.0 &&
      options.build.support_index_path.empty()) {
    throw std::runtime_error(
        "--minimum-supported-fraction requires --support-index");
  }
  return options;
}

VariantBitsMergeOptions ParseVariantBitsMergeOptions(int argc, char** argv) {
  if (argc == 3 && std::string(argv[2]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  VariantBitsMergeOptions options;
  std::unordered_map<std::string, std::string*> destinations{
      {"--list", &options.list_path},
      {"--out", &options.output_path},
  };
  for (int idx = 2; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    const auto iter = destinations.find(key);
    if (iter == destinations.end()) {
      throw std::runtime_error("unknown merge-variant-bits argument: " + key);
    }
    if (++idx >= argc) throw std::runtime_error("missing value after " + key);
    if (!iter->second->empty()) {
      throw std::runtime_error("argument supplied twice: " + key);
    }
    *iter->second = argv[idx];
  }
  if (options.list_path.empty() || options.output_path.empty()) {
    throw std::runtime_error("merge-variant-bits requires --list and --out");
  }
  return options;
}

SupportIndexOptions ParseSupportIndexOptions(int argc, char** argv) {
  if (argc == 3 && std::string(argv[2]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  SupportIndexOptions options;
  std::unordered_map<std::string, std::string*> destinations{
      {"--variant-index", &options.build.variant_index_path},
      {"--pvar", &options.build.pvar_path},
      {"--pvar-list", &options.build.pvar_list_path},
      {"--read-freq", &options.build.frequency_path},
      {"--progress-jsonl", &options.progress_jsonl},
      {"--progress-interval-seconds", &options.progress_interval_seconds},
      {"--out", &options.build.output_path},
  };
  for (int idx = 2; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    const auto iter = destinations.find(key);
    if (iter == destinations.end()) {
      throw std::runtime_error("unknown build-support-index argument: " + key);
    }
    if (++idx >= argc) throw std::runtime_error("missing value after " + key);
    if (!iter->second->empty()) {
      throw std::runtime_error("argument supplied twice: " + key);
    }
    *iter->second = argv[idx];
  }
  if (options.build.variant_index_path.empty() ||
      (options.build.pvar_path.empty() == options.build.pvar_list_path.empty()) ||
      options.build.frequency_path.empty() ||
      options.build.output_path.empty()) {
    throw std::runtime_error(
        "build-support-index requires --variant-index, exactly one of --pvar "
        "and --pvar-list, --read-freq, and --out");
  }
  if (!options.progress_interval_seconds.empty() &&
      options.progress_jsonl.empty()) {
    throw std::runtime_error(
        "--progress-interval-seconds requires --progress-jsonl");
  }
  if (!options.progress_interval_seconds.empty()) {
    ParseProgressInterval(options.progress_interval_seconds);
  }
  return options;
}

FragmentSupportOptions ParseFragmentSupportOptions(int argc, char** argv) {
  if (argc == 3 && std::string(argv[2]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  FragmentSupportOptions options;
  std::unordered_map<std::string, std::string*> destinations{
      {"--fragment", &options.fragment_path},
      {"--support-index", &options.support_index_path},
      {"--out", &options.output_path},
  };
  for (int idx = 2; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    const auto iter = destinations.find(key);
    if (iter == destinations.end()) {
      throw std::runtime_error("unknown report-fragment-support argument: " +
                               key);
    }
    if (++idx >= argc) throw std::runtime_error("missing value after " + key);
    if (!iter->second->empty()) {
      throw std::runtime_error("argument supplied twice: " + key);
    }
    *iter->second = argv[idx];
  }
  if (options.fragment_path.empty() || options.support_index_path.empty() ||
      options.output_path.empty()) {
    throw std::runtime_error(
        "report-fragment-support requires --fragment, --support-index, and "
        "--out");
  }
  return options;
}

std::unique_ptr<pgensparsescore::ProgressReporter> MakeProgressReporter(
    const std::string& path, const std::string& interval_seconds) {
  if (path.empty()) {
    return nullptr;
  }
  return std::make_unique<pgensparsescore::ProgressReporter>(
      path, ParseProgressInterval(interval_seconds));
}

pgensparsescore::MissingFrequencyPolicy ParseMissingFrequencyPolicy(
    const std::string& value) {
  if (value == "cohort") {
    return pgensparsescore::MissingFrequencyPolicy::kCohort;
  }
  if (value == "error") {
    return pgensparsescore::MissingFrequencyPolicy::kError;
  }
  if (value == "omit") {
    return pgensparsescore::MissingFrequencyPolicy::kOmit;
  }
  throw std::runtime_error(
      "--missing-freq must be cohort, error, or omit");
}

std::string MissingFrequencyPolicyName(
    pgensparsescore::MissingFrequencyPolicy policy) {
  switch (policy) {
    case pgensparsescore::MissingFrequencyPolicy::kCohort:
      return "cohort";
    case pgensparsescore::MissingFrequencyPolicy::kError:
      return "error";
    case pgensparsescore::MissingFrequencyPolicy::kOmit:
      return "omit";
  }
  throw std::runtime_error("invalid missing-frequency policy");
}

CompileOptions ParseCompileOptions(int argc, char** argv) {
  if (argc == 3 && std::string(argv[2]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  CompileOptions options;
  std::unordered_map<std::string, std::string*> destinations{
      {"--manifest", &options.manifest},
      {"--variant-map", &options.variant_map},
      {"--progress-jsonl", &options.progress_jsonl},
      {"--progress-interval-seconds", &options.progress_interval_seconds},
      {"--out", &options.out},
  };
  for (int idx = 2; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    const auto iter = destinations.find(key);
    if (iter == destinations.end()) {
      throw std::runtime_error("unknown compile argument: " + key);
    }
    if (++idx >= argc) {
      throw std::runtime_error("missing value after " + key);
    }
    if (!iter->second->empty()) {
      throw std::runtime_error("argument supplied twice: " + key);
    }
    *iter->second = argv[idx];
  }
  if (options.manifest.empty() || options.out.empty()) {
    throw std::runtime_error("compile requires --manifest and --out");
  }
  if (!options.progress_interval_seconds.empty() &&
      options.progress_jsonl.empty()) {
    throw std::runtime_error(
        "--progress-interval-seconds requires --progress-jsonl");
  }
  if (!options.progress_interval_seconds.empty()) {
    ParseProgressInterval(options.progress_interval_seconds);
  }
  if (options.progress_jsonl == options.out ||
      options.progress_jsonl == options.out + ".json") {
    throw std::runtime_error("progress log and compile output must differ");
  }
  return options;
}

Options ParseOptions(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  Options options;
  std::unordered_map<std::string, std::string*> destinations{
      {"--pgen", &options.pgen},
      {"--pvar", &options.pvar},
      {"--psam", &options.psam},
      {"--pfile-list", &options.pfile_list},
      {"--manifest", &options.manifest},
      {"--compiled-catalog", &options.compiled_catalog},
      {"--fragment-list", &options.fragment_list},
      {"--variant-index", &options.variant_index},
      {"--support-index", &options.support_index},
      {"--score-schema", &options.score_schema},
      {"--variant-map", &options.variant_map},
      {"--read-freq", &options.read_freq},
      {"--progress-jsonl", &options.progress_jsonl},
      {"--progress-interval-seconds", &options.progress_interval_seconds},
      {"--threads", &options.threads},
      {"--dense-kernel", &options.dense_kernel},
      {"--output-format", &options.output_format},
      {"--out", &options.out},
  };
  for (int idx = 1; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    if (key == "--error-on-missing-freq") {
      if (options.missing_frequency_policy_supplied) {
        throw std::runtime_error("argument supplied twice: " + key);
      }
      options.missing_frequency_policy =
          pgensparsescore::MissingFrequencyPolicy::kError;
      options.missing_frequency_policy_supplied = true;
      continue;
    }
    if (key == "--missing-freq") {
      if (options.missing_frequency_policy_supplied) {
        throw std::runtime_error("argument supplied twice: " + key);
      }
      if (++idx >= argc) {
        throw std::runtime_error("missing value after " + key);
      }
      options.missing_frequency_policy =
          ParseMissingFrequencyPolicy(argv[idx]);
      options.missing_frequency_policy_supplied = true;
      continue;
    }
    const auto iter = destinations.find(key);
    if (iter == destinations.end()) {
      throw std::runtime_error("unknown argument: " + key);
    }
    if (++idx >= argc) {
      throw std::runtime_error("missing value after " + key);
    }
    if (!iter->second->empty()) {
      throw std::runtime_error("argument supplied twice: " + key);
    }
    *iter->second = argv[idx];
  }
  const uint32_t score_source_ct = !options.manifest.empty() +
                                   !options.compiled_catalog.empty() +
                                   !options.fragment_list.empty();
  if (options.out.empty() || score_source_ct != 1) {
    throw std::runtime_error(
        "--out and exactly one of --manifest, --compiled-catalog, or "
        "--fragment-list are required");
  }
  if (!options.fragment_list.empty() && options.variant_index.empty()) {
    throw std::runtime_error("--fragment-list requires --variant-index");
  }
  if (options.fragment_list.empty() && !options.variant_index.empty()) {
    throw std::runtime_error("--variant-index requires --fragment-list");
  }
  if (options.fragment_list.empty() && !options.support_index.empty()) {
    throw std::runtime_error("--support-index requires --fragment-list");
  }
  if (options.fragment_list.empty() && !options.score_schema.empty()) {
    throw std::runtime_error("--score-schema requires --fragment-list");
  }
  if (options.fragment_list.empty() && !options.threads.empty()) {
    throw std::runtime_error("--threads requires --fragment-list");
  }
  if (options.fragment_list.empty() && !options.dense_kernel.empty()) {
    throw std::runtime_error("--dense-kernel requires --fragment-list");
  }
  if (!options.fragment_list.empty() && !options.variant_map.empty()) {
    throw std::runtime_error(
        "fragment indexes already contain ID aliases; do not use --variant-map");
  }
  const bool has_single = !options.pgen.empty() || !options.pvar.empty() ||
                          !options.psam.empty();
  if (!options.pfile_list.empty() && has_single) {
    throw std::runtime_error(
        "--pfile-list cannot be combined with --pgen/--pvar/--psam");
  }
  if (options.pfile_list.empty() &&
      (options.pgen.empty() || options.pvar.empty() || options.psam.empty())) {
    throw std::runtime_error(
        "supply either --pfile-list or all of --pgen/--pvar/--psam");
  }
  if (options.missing_frequency_policy !=
          pgensparsescore::MissingFrequencyPolicy::kCohort &&
      options.read_freq.empty()) {
    throw std::runtime_error(
        "error and omit missing-frequency policies require --read-freq");
  }
  if (!options.progress_interval_seconds.empty() &&
      options.progress_jsonl.empty()) {
    throw std::runtime_error(
        "--progress-interval-seconds requires --progress-jsonl");
  }
  if (!options.progress_interval_seconds.empty()) {
    ParseProgressInterval(options.progress_interval_seconds);
  }
  if (!options.threads.empty()) {
    ParsePositiveU32(options.threads, "--threads");
  }
  const auto dense_kernel = ParseDenseScoringKernel(options.dense_kernel);
  if (dense_kernel == pgensparsescore::DenseScoringKernel::kOneMkl &&
      !pgensparsescore::OneMklDenseScoringAvailable()) {
    throw std::runtime_error(
        "--dense-kernel onemkl requires a oneMKL-enabled build");
  }
  ParseScoreOutputFormat(options.output_format);
  return options;
}

void AppendDouble(double value, std::string* output) {
  char buffer[64];
  const auto result = std::to_chars(
      buffer, buffer + sizeof(buffer), value, std::chars_format::general,
      std::numeric_limits<double>::max_digits10);
  if (result.ec != std::errc()) {
    throw std::runtime_error("cannot format score value");
  }
  output->append(buffer, result.ptr);
}

void WriteWideScores(const std::string& prefix,
                     const std::vector<pgensparsescore::Sample>& samples,
                     const pgensparsescore::Catalog& catalog,
                     const pgensparsescore::MappedMatrix& matrix,
                     pgensparsescore::ProgressReporter* progress) {
  const std::string final_path = prefix + ".scores.tsv.gz";
  const std::string temporary_path = final_path + ".tmp";
  RemoveFileOnExit remove_temporary(temporary_path);
  pgensparsescore::GzipWriter output(temporary_path);

  const bool has_fid = samples.front().fid.has_value();
  std::string line = has_fid ? "FID\tIID" : "IID";
  for (const auto& score : catalog.scores) {
    line.push_back('\t');
    line += score.id;
  }
  line.push_back('\n');
  output.Write(line);

  constexpr uint64_t kTransposeBufferBytes = 64ULL * 1024 * 1024;
  const uint64_t bytes_per_sample =
      static_cast<uint64_t>(catalog.scores.size()) * sizeof(double);
  const uint32_t block_capacity = static_cast<uint32_t>(std::max<uint64_t>(
      1, std::min<uint64_t>(samples.size(),
                            kTransposeBufferBytes / bytes_per_sample)));
  std::vector<double> sample_major_block(
      static_cast<uint64_t>(block_capacity) * catalog.scores.size());

  if (progress) {
    progress->Event("score", "write_scores",
                    {{"sample_rows_written", 0},
                     {"sample_rows_total", samples.size()},
                     {"score_columns", catalog.scores.size()}});
  }

  for (uint32_t sample_begin = 0; sample_begin < samples.size();
       sample_begin += block_capacity) {
    const uint32_t block_size = static_cast<uint32_t>(std::min<uint64_t>(
        block_capacity, samples.size() - sample_begin));
    for (uint32_t score_idx = 0; score_idx < catalog.scores.size();
         ++score_idx) {
      const double* score_row = matrix.Row(score_idx) + sample_begin;
      for (uint32_t block_idx = 0; block_idx < block_size; ++block_idx) {
        sample_major_block[static_cast<uint64_t>(block_idx) *
                               catalog.scores.size() +
                           score_idx] = score_row[block_idx];
      }
    }
    for (uint32_t block_idx = 0; block_idx < block_size; ++block_idx) {
      const auto& sample = samples[sample_begin + block_idx];
      line.clear();
      line.reserve(32 + catalog.scores.size() * 12);
      if (has_fid) {
        line += *sample.fid;
        line.push_back('\t');
      }
      line += sample.iid;
      const double* row =
          sample_major_block.data() +
          static_cast<uint64_t>(block_idx) * catalog.scores.size();
      for (uint32_t score_idx = 0; score_idx < catalog.scores.size();
           ++score_idx) {
        line.push_back('\t');
        AppendDouble(row[score_idx], &line);
      }
      line.push_back('\n');
      output.Write(line);
    }
    if (progress) {
      progress->MaybeEvent(
          "score", "write_scores",
          {{"sample_rows_written", sample_begin + block_size},
           {"sample_rows_total", samples.size()},
           {"score_columns", catalog.scores.size()}});
    }
  }
  output.Close();
  std::filesystem::rename(temporary_path, final_path);
  remove_temporary.Release();
  if (progress) {
    progress->Event("score", "scores_written",
                    {{"sample_rows_written", samples.size()},
                     {"score_columns", catalog.scores.size()},
                     {"output_bytes", std::filesystem::file_size(final_path)}});
  }
}

void WriteSampleTable(const std::string& prefix,
                      const std::vector<pgensparsescore::Sample>& samples) {
  const std::string final_path = prefix + ".samples.tsv";
  const std::string temporary_path = final_path + ".tmp";
  RemoveFileOnExit remove_temporary(temporary_path);
  std::ofstream output(temporary_path);
  if (!output) throw std::runtime_error("cannot write sample table");
  const bool has_fid = samples.front().fid.has_value();
  output << (has_fid ? "FID\tIID\n" : "IID\n");
  for (const auto& sample : samples) {
    if (has_fid) output << *sample.fid << '\t';
    output << sample.iid << '\n';
  }
  output.close();
  if (!output) throw std::runtime_error("cannot finish sample table");
  std::filesystem::rename(temporary_path, final_path);
  remove_temporary.Release();
}

void FinishBinaryScores(
    const std::string& prefix,
    const std::vector<pgensparsescore::Sample>& samples,
    const pgensparsescore::Catalog& catalog,
    pgensparsescore::MappedMatrix* matrix, const std::string& working_path,
    RemoveFileOnExit* remove_working,
    pgensparsescore::ProgressReporter* progress) {
  const uint16_t endian_probe = 1;
  if (*reinterpret_cast<const uint8_t*>(&endian_probe) != 1) {
    throw std::runtime_error(
        "score-major binary output requires a little-endian machine");
  }
  matrix->Flush();
  WriteSampleTable(prefix, samples);
  const std::string final_path = ScoreOutputPath(
      prefix, ScoreOutputFormat::kScoreMajorBinary);
  std::filesystem::rename(working_path, final_path);
  remove_working->Release();
  if (progress) {
    progress->Event("score", "scores_written",
                    {{"sample_rows_written", samples.size()},
                     {"score_columns", catalog.scores.size()},
                     {"output_bytes", std::filesystem::file_size(final_path)}},
                    {{"output_format", "score-major-bin"}});
  }
}

bool SamplesEqual(const std::vector<pgensparsescore::Sample>& lhs,
                  const std::vector<pgensparsescore::Sample>& rhs) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (size_t idx = 0; idx < lhs.size(); ++idx) {
    if (lhs[idx].iid != rhs[idx].iid || lhs[idx].fid != rhs[idx].fid) {
      return false;
    }
  }
  return true;
}

std::vector<pgensparsescore::Catalog> PartitionCatalog(
    pgensparsescore::Catalog* catalog,
    const std::vector<uint32_t>& input_by_variant,
    const std::vector<uint32_t>& local_index_by_variant,
    size_t input_ct) {
  std::vector<pgensparsescore::Catalog> result(input_ct);
  for (auto& shard : result) {
    shard.scores.resize(catalog->scores.size());
    shard.intercepts.assign(catalog->scores.size(), 0.0);
  }
  result.front().intercepts = catalog->intercepts;
  for (auto& variant : catalog->variants) {
    const uint32_t input_idx = input_by_variant.at(variant.variant_idx);
    result.at(input_idx).variants.push_back(
        {local_index_by_variant.at(variant.variant_idx),
         std::move(variant.edges)});
  }
  catalog->variants.clear();
  catalog->variants.shrink_to_fit();
  return result;
}

void AddStats(const pgensparsescore::ScoreRunStats& input,
              pgensparsescore::ScoreRunStats* output) {
  output->variant_ct += input.variant_ct;
  output->edge_ct += input.edge_ct;
  output->sparse_variant_ct += input.sparse_variant_ct;
  output->dense_variant_ct += input.dense_variant_ct;
  output->sparse_edge_ct += input.sparse_edge_ct;
  output->dense_edge_ct += input.dense_edge_ct;
  output->sparse_value_ct += input.sparse_value_ct;
  output->sparse_update_ct += input.sparse_update_ct;
  output->dense_update_ct += input.dense_update_ct;
  output->parallel_variant_ct += input.parallel_variant_ct;
  output->parallel_update_ct += input.parallel_update_ct;
  output->score_major_tile_ct += input.score_major_tile_ct;
  output->score_major_row_ct += input.score_major_row_ct;
  output->score_major_maximum_rows_per_tile =
      std::max(output->score_major_maximum_rows_per_tile,
               input.score_major_maximum_rows_per_tile);
  output->score_major_maximum_edges_per_tile =
      std::max(output->score_major_maximum_edges_per_tile,
               input.score_major_maximum_edges_per_tile);
  output->score_major_scoring_nanoseconds +=
      input.score_major_scoring_nanoseconds;
  output->direct_dense_tile_ct += input.direct_dense_tile_ct;
  output->onemkl_tile_ct += input.onemkl_tile_ct;
  output->onemkl_matrix_build_nanoseconds +=
      input.onemkl_matrix_build_nanoseconds;
  output->onemkl_optimize_nanoseconds += input.onemkl_optimize_nanoseconds;
  output->onemkl_multiply_nanoseconds += input.onemkl_multiply_nanoseconds;
  output->densified_sparse_variant_ct += input.densified_sparse_variant_ct;
  output->copied_sparse_genotype_bytes += input.copied_sparse_genotype_bytes;
  output->maximum_genotype_buffer_bytes =
      std::max(output->maximum_genotype_buffer_bytes,
               input.maximum_genotype_buffer_bytes);
  output->imputed_value_ct += input.imputed_value_ct;
  output->external_frequency_variant_ct += input.external_frequency_variant_ct;
  output->cohort_frequency_variant_ct += input.cohort_frequency_variant_ct;
  output->missing_frequency_variant_ct += input.missing_frequency_variant_ct;
  output->omitted_frequency_variant_ct += input.omitted_frequency_variant_ct;
  output->omitted_frequency_edge_ct += input.omitted_frequency_edge_ct;
}

void WriteMetadata(const std::string& prefix, uint32_t sample_ct, bool has_fid,
                   uint32_t pgen_ct, uint64_t frequency_row_ct,
                   uint64_t variant_mapping_row_ct,
                   uint64_t pvar_variant_ct,
                   const std::vector<std::string>& pgen_storage_modes,
                   const std::string& catalog_source,
                   pgensparsescore::MissingFrequencyPolicy frequency_policy,
                   uint64_t working_matrix_byte_ct,
                   const pgensparsescore::Catalog& catalog,
                   const pgensparsescore::ScoreRunStats& stats,
                   uint32_t score_fragment_ct = 0,
                   uint64_t variant_index_variant_ct = 0,
                   uint32_t scoring_thread_ct = 1,
                   pgensparsescore::DenseScoringKernel dense_kernel =
                       pgensparsescore::DenseScoringKernel::kDirect,
                   uint32_t onemkl_thread_ct = 0,
                   ScoreOutputFormat output_format =
                       ScoreOutputFormat::kWideTsv) {
  {
    std::ofstream output(prefix + ".score-metadata.tsv");
    if (!output) throw std::runtime_error("cannot write score metadata");
    output << "INDEX\tSCORE\tINPUT_WEIGHTS\tZERO_WEIGHTS"
              "\tEXCLUDED_WEIGHTS\tDUPLICATE_WEIGHTS\tCATALOG_WEIGHTS\tMATCHED_WEIGHTS"
              "\tMISSING_VARIANTS\tMISSING_FREQUENCIES\tSCORED_WEIGHTS"
              "\tALT_EFFECTS\tREF_EFFECTS\tNONZERO_WEIGHT_L1"
              "\tNONZERO_WEIGHT_L2_SQUARED\tCATALOG_WEIGHT_L1"
              "\tCATALOG_WEIGHT_L2_SQUARED"
              "\tREFERENCE_SUPPORTED_WEIGHT_L1"
              "\tREFERENCE_SUPPORTED_WEIGHT_L2_SQUARED\tREF_INTERCEPT\n";
    for (uint32_t idx = 0; idx < catalog.scores.size(); ++idx) {
      const auto& score = catalog.scores[idx];
      output << idx << '\t' << score.id << '\t' << score.input_weight_ct << '\t'
             << score.zero_weight_ct << '\t' << score.excluded_weight_ct << '\t'
             << score.duplicate_weight_ct << '\t' << score.catalog_weight_ct
             << '\t' << score.matched_weight_ct
             << '\t' << score.missing_variant_ct << '\t'
             << score.missing_frequency_ct << '\t'
             << (score.matched_weight_ct - score.missing_frequency_ct)
             << '\t' << score.alt_effect_ct << '\t' << score.ref_effect_ct
             << '\t' << score.nonzero_weight_l1 << '\t'
             << score.nonzero_weight_l2 << '\t' << score.catalog_weight_l1
             << '\t' << score.catalog_weight_l2 << '\t'
             << score.supported_weight_l1 << '\t'
             << score.supported_weight_l2 << '\t'
             << score.ref_effect_intercept << '\n';
    }
  }
  {
    std::ofstream output(prefix + ".json");
    if (!output) throw std::runtime_error("cannot write JSON metadata");
    const std::filesystem::path scores(
        ScoreOutputPath(prefix, output_format));
    const char* dense_kernel_used = "direct";
    if (score_fragment_ct) {
      if (stats.direct_dense_tile_ct && stats.onemkl_tile_ct) {
        dense_kernel_used = "mixed";
      } else if (stats.onemkl_tile_ct) {
        dense_kernel_used = "onemkl";
      } else if (!stats.direct_dense_tile_ct) {
        dense_kernel_used = "none";
      }
    }
    output << "{\n"
           << "  \"format\": \"" << ScoreOutputFormatName(output_format)
           << "\",\n"
           << "  \"path\": \""
           << pgensparsescore::JsonEscape(scores.filename().string()) << "\",\n";
    if (output_format == ScoreOutputFormat::kScoreMajorBinary) {
      output << "  \"samples_path\": \""
             << pgensparsescore::JsonEscape(
                    std::filesystem::path(prefix + ".samples.tsv")
                        .filename()
                        .string())
             << "\",\n"
                "  \"dtype\": \"float64\",\n"
                "  \"byte_order\": \"little-endian\",\n"
                "  \"matrix_layout\": \"score-major\",\n";
    }
    output << "  \"sample_id_columns\": "
           << (has_fid ? "[\"FID\", \"IID\"]" : "[\"IID\"]") << ",\n"
           << "  \"sample_rows\": " << sample_ct << ",\n"
           << "  \"score_columns\": " << catalog.scores.size() << ",\n"
           << "  \"pgen_inputs\": " << pgen_ct << ",\n"
           << "  \"frequency_rows\": " << frequency_row_ct << ",\n"
           << "  \"variant_mapping_rows\": " << variant_mapping_row_ct
           << ",\n"
           << "  \"pvar_variants_loaded\": " << pvar_variant_ct << ",\n"
           << "  \"pgen_storage_modes\": [";
    for (size_t idx = 0; idx < pgen_storage_modes.size(); ++idx) {
      if (idx) output << ", ";
      output << "\"" << pgen_storage_modes[idx] << "\"";
    }
    output << "],\n"
           << "  \"catalog_source\": \"" << catalog_source << "\",\n"
           << "  \"score_fragments\": " << score_fragment_ct << ",\n"
           << "  \"variant_index_variants\": "
           << variant_index_variant_ct << ",\n"
           << "  \"missing_frequency_policy\": \""
           << MissingFrequencyPolicyName(frequency_policy) << "\",\n"
           << "  \"working_matrix_bytes\": " << working_matrix_byte_ct
           << ",\n"
           << "  \"scoring_threads\": " << scoring_thread_ct << ",\n"
           << "  \"onemkl_threads\": " << onemkl_thread_ct << ",\n"
           << "  \"dense_scoring_kernel_requested\": \""
           << pgensparsescore::DenseScoringKernelName(dense_kernel)
           << "\",\n"
           << "  \"dense_scoring_kernel_used\": \""
           << dense_kernel_used << "\",\n"
           << "  \"scoring_layout\": \""
           << (score_fragment_ct ? "score-major-tiles" : "variant-major")
           << "\",\n"
           << "  \"scored_variants\": " << stats.variant_ct << ",\n"
           << "  \"weight_edges\": " << stats.edge_ct << ",\n"
           << "  \"sparse_variants\": " << stats.sparse_variant_ct << ",\n"
           << "  \"dense_variants\": " << stats.dense_variant_ct << ",\n"
           << "  \"sparse_weight_edges\": " << stats.sparse_edge_ct << ",\n"
           << "  \"dense_weight_edges\": " << stats.dense_edge_ct << ",\n"
           << "  \"sparse_dosage_values\": " << stats.sparse_value_ct << ",\n"
           << "  \"sparse_score_updates\": " << stats.sparse_update_ct << ",\n"
           << "  \"dense_score_updates\": " << stats.dense_update_ct << ",\n"
           << "  \"parallel_variants\": " << stats.parallel_variant_ct
           << ",\n"
           << "  \"parallel_score_updates\": " << stats.parallel_update_ct
           << ",\n"
           << "  \"score_major_tiles\": " << stats.score_major_tile_ct
           << ",\n"
           << "  \"score_major_rows\": " << stats.score_major_row_ct
           << ",\n"
           << "  \"score_major_maximum_rows_per_tile\": "
           << stats.score_major_maximum_rows_per_tile << ",\n"
           << "  \"score_major_maximum_edges_per_tile\": "
           << stats.score_major_maximum_edges_per_tile << ",\n"
           << "  \"score_major_scoring_nanoseconds\": "
           << stats.score_major_scoring_nanoseconds << ",\n"
           << "  \"direct_dense_tiles\": " << stats.direct_dense_tile_ct
           << ",\n"
           << "  \"onemkl_tiles\": " << stats.onemkl_tile_ct << ",\n"
           << "  \"onemkl_matrix_build_nanoseconds\": "
           << stats.onemkl_matrix_build_nanoseconds << ",\n"
           << "  \"onemkl_optimize_nanoseconds\": "
           << stats.onemkl_optimize_nanoseconds << ",\n"
           << "  \"onemkl_multiply_nanoseconds\": "
           << stats.onemkl_multiply_nanoseconds << ",\n"
           << "  \"densified_sparse_variants\": "
           << stats.densified_sparse_variant_ct << ",\n"
           << "  \"copied_sparse_genotype_bytes\": "
           << stats.copied_sparse_genotype_bytes << ",\n"
           << "  \"maximum_genotype_buffer_bytes\": "
           << stats.maximum_genotype_buffer_bytes << ",\n"
           << "  \"imputed_values\": " << stats.imputed_value_ct << ",\n"
           << "  \"external_frequency_variants\": "
           << stats.external_frequency_variant_ct << ",\n"
           << "  \"cohort_frequency_variants\": "
           << stats.cohort_frequency_variant_ct << ",\n"
           << "  \"missing_frequency_variants\": "
           << stats.missing_frequency_variant_ct << ",\n"
           << "  \"omitted_frequency_variants\": "
           << stats.omitted_frequency_variant_ct << ",\n"
           << "  \"omitted_frequency_edges\": "
           << stats.omitted_frequency_edge_ct << "\n"
           << "}\n";
  }
}

void WriteCompiledCatalogMetadata(
    const std::string& path,
    const pgensparsescore::CompiledCatalog& catalog,
    uint64_t included_source_variant_ct) {
  std::ofstream output(path + ".json");
  if (!output) {
    throw std::runtime_error("cannot write compiled-catalog metadata");
  }
  uint64_t duplicate_weight_ct = 0;
  for (const auto& score : catalog.scores) {
    duplicate_weight_ct += score.duplicate_weight_ct;
  }
  output << "{\n"
         << "  \"format\": \"pgensparsescore-compiled-catalog-v2\",\n"
         << "  \"path\": \""
         << pgensparsescore::JsonEscape(
                std::filesystem::path(path).filename().string())
         << "\",\n"
         << "  \"scores\": " << catalog.scores.size() << ",\n"
         << "  \"variants\": " << catalog.variants.size() << ",\n"
         << "  \"weights\": " << catalog.weight_ct << ",\n"
         << "  \"duplicate_weights\": " << duplicate_weight_ct << ",\n"
         << "  \"included_source_variants\": "
         << included_source_variant_ct << "\n"
         << "}\n";
}

int RunFragmentScoring(
    const Options& options,
    const std::vector<pgensparsescore::PfileSpec>& inputs,
    pgensparsescore::ProgressReporter* progress) {
  pgensparsescore::VariantIndex variant_index(options.variant_index);
  std::unique_ptr<pgensparsescore::SupportIndex> support_index;
  if (!options.support_index.empty()) {
    support_index =
        std::make_unique<pgensparsescore::SupportIndex>(options.support_index);
    if (support_index->variant_ct() != variant_index.variant_ct() ||
        support_index->signature_lo() != variant_index.signature_lo() ||
        support_index->signature_hi() != variant_index.signature_hi()) {
      throw std::runtime_error(
          "support index was built for a different variant index");
    }
  }
  const auto fragment_paths =
      pgensparsescore::ReadScoreFragmentList(options.fragment_list);
  auto loaded = pgensparsescore::LoadScoreFragments(
      fragment_paths, variant_index, options.score_schema, progress);
  auto samples = pgensparsescore::ReadPsam(inputs.front().psam);
  if (progress) {
    progress->Event("score", "samples_loaded", {{"samples", samples.size()}});
  }

  std::vector<pgensparsescore::IndexedVariantLocation> locations(
      variant_index.variant_ct());
  std::vector<uint32_t> pvar_row_counts;
  pvar_row_counts.reserve(inputs.size());
  uint64_t indexed_pvar_variant_ct = 0;
  for (uint32_t input_idx = 0; input_idx < inputs.size(); ++input_idx) {
    if (input_idx) {
      const auto input_samples =
          pgensparsescore::ReadPsam(inputs[input_idx].psam);
      if (!SamplesEqual(samples, input_samples)) {
        throw std::runtime_error(
            "PSAM sample IDs/order differ between PGEN inputs: " +
            inputs[input_idx].psam);
      }
    }
    const auto pvar = pgensparsescore::AddIndexedPvar(
        inputs[input_idx].pvar, input_idx, variant_index, &locations, progress);
    pvar_row_counts.push_back(pvar.row_ct);
    indexed_pvar_variant_ct += pvar.matched_variant_ct;
  }

  std::optional<pgensparsescore::IndexedFrequencyTable> frequencies;
  if (!options.read_freq.empty()) {
    frequencies = pgensparsescore::ReadIndexedFrequencyTable(
        options.read_freq, variant_index, progress);
  }

  std::vector<std::unique_ptr<pgensparsescore::PgenDosageReader>> owned_readers;
  std::vector<pgensparsescore::PgenDosageReader*> readers;
  std::vector<std::string> storage_modes;
  owned_readers.reserve(inputs.size());
  readers.reserve(inputs.size());
  storage_modes.reserve(inputs.size());
  for (uint32_t input_idx = 0; input_idx < inputs.size(); ++input_idx) {
    owned_readers.push_back(
        std::make_unique<pgensparsescore::PgenDosageReader>(
            inputs[input_idx].pgen));
    auto* reader = owned_readers.back().get();
    if (reader->variant_ct() != pvar_row_counts[input_idx]) {
      throw std::runtime_error("PGEN/PVAR variant-count mismatch: " +
                               inputs[input_idx].pgen);
    }
    if (reader->sample_ct() != samples.size()) {
      throw std::runtime_error("PGEN/PSAM sample-count mismatch: " +
                               inputs[input_idx].pgen);
    }
    readers.push_back(reader);
    storage_modes.emplace_back(reader->storage_mode_name());
    if (progress) {
      progress->Event(
          "score", "pgen_ready",
          {{"pgen_input_index", input_idx},
           {"pgen_inputs", inputs.size()},
           {"pgen_variants", reader->variant_ct()},
           {"samples", reader->sample_ct()}},
          {{"pgen", inputs[input_idx].pgen},
           {"storage_mode", reader->storage_mode_name()}});
    }
  }

  const auto output_format = ParseScoreOutputFormat(options.output_format);
  const std::string working_path =
      output_format == ScoreOutputFormat::kScoreMajorBinary
          ? ScoreOutputPath(options.out, output_format) + ".tmp"
          : options.out + ".work.score-major.bin";
  RemoveFileOnExit remove_working(working_path);
  uint64_t working_matrix_byte_ct = 0;
  const uint32_t scoring_thread_ct = options.threads.empty()
                                         ? LogicalCoreCount()
                                         : ParsePositiveU32(options.threads,
                                                            "--threads");
  const auto dense_kernel = ParseDenseScoringKernel(options.dense_kernel);
  const uint32_t onemkl_thread_ct =
      dense_kernel != pgensparsescore::DenseScoringKernel::kDirect &&
              pgensparsescore::OneMklDenseScoringAvailable()
          ? (options.threads.empty() ? LogicalCoreCount() : scoring_thread_ct)
          : 0;
  pgensparsescore::ScoreRunStats stats;
  {
    pgensparsescore::MappedMatrix matrix(
        working_path, loaded.catalog.scores.size(), samples.size());
    working_matrix_byte_ct = matrix.byte_ct();
    if (progress) {
      progress->Event(
          "score", "working_matrix_ready",
          {{"working_matrix_bytes", working_matrix_byte_ct},
           {"sample_rows", samples.size()},
           {"score_columns", loaded.catalog.scores.size()},
           {"score_fragments", loaded.fragments.size()},
           {"fragment_weights", loaded.weight_ct},
           {"scoring_threads", scoring_thread_ct},
           {"onemkl_threads", onemkl_thread_ct}},
          {{"dense_scoring_kernel_requested",
            pgensparsescore::DenseScoringKernelName(dense_kernel)}});
    }
    stats = pgensparsescore::ScoreFragments(
        variant_index, loaded.fragments, loaded.score_maps, locations,
        frequencies ? &*frequencies : nullptr,
        options.missing_frequency_policy, readers, &loaded.catalog, &matrix,
        scoring_thread_ct, dense_kernel, onemkl_thread_ct, progress);
    if (support_index) {
      stats.missing_frequency_variant_ct +=
          support_index->missing_frequency_ct();
      stats.omitted_frequency_variant_ct +=
          support_index->missing_frequency_ct();
    }
    if (output_format == ScoreOutputFormat::kWideTsv) {
      WriteWideScores(options.out, samples, loaded.catalog, matrix, progress);
    } else {
      FinishBinaryScores(options.out, samples, loaded.catalog, &matrix,
                         working_path, &remove_working, progress);
    }
  }
  if (output_format == ScoreOutputFormat::kWideTsv) {
    remove_working.RemoveNow();
  }
  WriteMetadata(
      options.out, samples.size(), samples.front().fid.has_value(),
      inputs.size(), frequencies ? frequencies->matched_row_ct : 0,
      variant_index.variant_ct(), indexed_pvar_variant_ct, storage_modes,
      "fragments", options.missing_frequency_policy, working_matrix_byte_ct,
      loaded.catalog, stats, loaded.fragments.size(), variant_index.variant_ct(),
      scoring_thread_ct, dense_kernel, onemkl_thread_ct, output_format);
  if (progress) {
    progress->Event(
        "score", "complete",
        {{"sample_rows", samples.size()},
         {"score_columns", loaded.catalog.scores.size()},
         {"score_fragments", loaded.fragments.size()},
         {"pgen_inputs", inputs.size()},
         {"genotype_decodes", stats.variant_ct},
         {"weight_edges", stats.edge_ct},
         {"sparse_variants", stats.sparse_variant_ct},
         {"dense_variants", stats.dense_variant_ct},
         {"sparse_weight_edges", stats.sparse_edge_ct},
         {"dense_weight_edges", stats.dense_edge_ct},
         {"sparse_score_updates", stats.sparse_update_ct},
         {"dense_score_updates", stats.dense_update_ct},
         {"score_major_tiles", stats.score_major_tile_ct},
         {"score_major_rows", stats.score_major_row_ct},
         {"score_major_scoring_nanoseconds",
          stats.score_major_scoring_nanoseconds},
         {"direct_dense_tiles", stats.direct_dense_tile_ct},
         {"densified_sparse_variants",
          stats.densified_sparse_variant_ct},
         {"maximum_genotype_buffer_bytes",
          stats.maximum_genotype_buffer_bytes},
         {"onemkl_tiles", stats.onemkl_tile_ct},
         {"onemkl_matrix_build_nanoseconds",
          stats.onemkl_matrix_build_nanoseconds},
         {"onemkl_optimize_nanoseconds",
          stats.onemkl_optimize_nanoseconds},
         {"onemkl_multiply_nanoseconds",
          stats.onemkl_multiply_nanoseconds},
         {"scoring_threads", scoring_thread_ct},
         {"onemkl_threads", onemkl_thread_ct},
        {"imputed_values", stats.imputed_value_ct}},
        {{"dense_scoring_kernel_requested",
          pgensparsescore::DenseScoringKernelName(dense_kernel)},
         {"output_format", ScoreOutputFormatName(output_format)}});
  }
  std::cerr << "wrote " << loaded.catalog.scores.size()
            << " named score columns for " << samples.size() << " samples from "
            << loaded.fragments.size() << " score fragments with "
            << stats.variant_ct << " total genotype decodes\n";
  return 0;
}

int CompileMain(int argc, char** argv) {
  const CompileOptions options = ParseCompileOptions(argc, argv);
  if (std::filesystem::exists(options.out) ||
      std::filesystem::exists(options.out + ".json")) {
    throw std::runtime_error("compile output already exists: " + options.out);
  }
  auto progress = MakeProgressReporter(options.progress_jsonl,
                                       options.progress_interval_seconds);
  if (progress) {
    progress->Event("compile", "start", {},
                    {{"manifest", options.manifest}, {"output", options.out}});
  }
  std::optional<pgensparsescore::VariantMap> variant_map;
  std::unordered_set<std::string> included_source_ids;
  if (!options.variant_map.empty()) {
    variant_map = pgensparsescore::ReadVariantMap(options.variant_map);
    included_source_ids.reserve(variant_map->size());
    for (const auto& [source_id, target_id] : *variant_map) {
      static_cast<void>(target_id);
      included_source_ids.insert(source_id);
    }
    if (progress) {
      progress->Event("compile", "variant_filter_loaded",
                      {{"included_variants", included_source_ids.size()}});
    }
  }
  const auto catalog = pgensparsescore::CompileSourceCatalog(
      options.manifest,
      variant_map ? &included_source_ids : nullptr, progress.get());
  pgensparsescore::WriteCompiledCatalog(options.out, catalog, progress.get());
  WriteCompiledCatalogMetadata(
      options.out, catalog,
      variant_map ? variant_map->size() : catalog.variants.size());
  if (progress) {
    progress->Event("compile", "complete",
                    {{"scores", catalog.scores.size()},
                     {"variants", catalog.variants.size()},
                     {"weights", catalog.weight_ct},
                     {"output_bytes", std::filesystem::file_size(options.out)}});
  }
  std::cerr << "compiled " << catalog.scores.size() << " scores, "
            << catalog.variants.size() << " variants, and "
            << catalog.weight_ct << " nonzero weights into " << options.out
            << '\n';
  return 0;
}

int VariantIndexMain(int argc, char** argv) {
  const VariantIndexOptions options = ParseVariantIndexOptions(argc, argv);
  if (std::filesystem::exists(options.build.output_path) ||
      std::filesystem::exists(options.build.output_path + ".json")) {
    throw std::runtime_error("variant-index output already exists: " +
                             options.build.output_path);
  }
  auto progress = MakeProgressReporter(options.progress_jsonl,
                                       options.progress_interval_seconds);
  if (progress) {
    progress->Event("variant_index", "start", {},
                    {{"variant_list", options.build.input_path},
                     {"output", options.build.output_path}});
  }
  pgensparsescore::BuildVariantIndex(options.build, progress.get());
  pgensparsescore::VariantIndex index(options.build.output_path);
  std::ofstream metadata(options.build.output_path + ".json");
  if (!metadata) {
    throw std::runtime_error("cannot write variant-index metadata");
  }
  metadata << "{\n"
           << "  \"format\": \"pgensparsescore-variant-index-v1\",\n"
           << "  \"path\": \""
           << pgensparsescore::JsonEscape(
                  std::filesystem::path(options.build.output_path)
                      .filename()
                      .string())
           << "\",\n"
           << "  \"variants\": " << index.variant_ct() << ",\n"
           << "  \"aliases\": " << index.alias_ct() << ",\n"
           << "  \"block_size\": " << index.block_size() << ",\n"
           << "  \"blocks\": " << index.block_ct() << ",\n"
           << "  \"bytes\": "
           << std::filesystem::file_size(options.build.output_path) << "\n"
           << "}\n";
  metadata.close();
  if (!metadata) throw std::runtime_error("cannot finish variant-index metadata");
  std::cerr << "built index for " << index.variant_ct() << " variants and "
            << index.alias_ct() << " IDs in " << index.block_ct()
            << " blocks\n";
  return 0;
}

int SupportIndexMain(int argc, char** argv) {
  const SupportIndexOptions options = ParseSupportIndexOptions(argc, argv);
  if (std::filesystem::exists(options.build.output_path) ||
      std::filesystem::exists(options.build.output_path + ".json")) {
    throw std::runtime_error("support-index output already exists: " +
                             options.build.output_path);
  }
  auto progress = MakeProgressReporter(options.progress_jsonl,
                                       options.progress_interval_seconds);
  const auto summary =
      pgensparsescore::BuildSupportIndex(options.build, progress.get());
  std::ofstream metadata(options.build.output_path + ".json");
  if (!metadata) throw std::runtime_error("cannot write support-index metadata");
  metadata << "{\n"
           << "  \"format\": \"pgensparsescore-support-index-v1\",\n"
           << "  \"path\": \""
           << pgensparsescore::JsonEscape(
                  std::filesystem::path(options.build.output_path)
                      .filename()
                      .string())
           << "\",\n"
           << "  \"variants\": " << summary.variant_ct << ",\n"
           << "  \"missing_variants\": " << summary.missing_variant_ct
           << ",\n"
           << "  \"missing_frequencies\": "
           << summary.missing_frequency_ct << ",\n"
           << "  \"usable_variants\": " << summary.usable_variant_ct
           << ",\n"
           << "  \"pvar_inputs\": " << summary.pvar_input_ct << ",\n"
           << "  \"pvar_rows\": " << summary.pvar_row_ct << ",\n"
           << "  \"frequency_rows\": " << summary.frequency_row_ct
           << ",\n"
           << "  \"bytes\": " << summary.output_bytes << "\n"
           << "}\n";
  metadata.close();
  if (!metadata) throw std::runtime_error("cannot finish support-index metadata");
  std::cerr << "built support index with " << summary.usable_variant_ct
            << " usable variants, " << summary.missing_variant_ct
            << " missing variants, and " << summary.missing_frequency_ct
            << " missing frequencies\n";
  return 0;
}

int FragmentSupportMain(int argc, char** argv) {
  const FragmentSupportOptions options =
      ParseFragmentSupportOptions(argc, argv);
  const std::string metadata_path = options.output_path + ".json";
  if (std::filesystem::exists(options.output_path) ||
      std::filesystem::exists(metadata_path)) {
    throw std::runtime_error("fragment-support output already exists: " +
                             options.output_path);
  }
  const pgensparsescore::ScoreFragmentReader fragment(options.fragment_path);
  const pgensparsescore::SupportIndex support(options.support_index_path);
  const auto summary =
      pgensparsescore::MeasureScoreFragmentSupport(fragment, support);

  const std::string temporary_output = options.output_path + ".tmp";
  RemoveFileOnExit output_cleanup(temporary_output);
  std::ofstream output(temporary_output);
  if (!output) throw std::runtime_error("cannot write " + temporary_output);
  output << "SCORE_ID\tCOLUMN_NAME\tREFERENCE_WEIGHT_COUNT"
            "\tTARGET_AVAILABLE_WEIGHT_COUNT\tTARGET_AVAILABLE_FRACTION"
            "\tREFERENCE_ABS_WEIGHT\tTARGET_AVAILABLE_ABS_WEIGHT"
            "\tTARGET_AVAILABLE_ABS_WEIGHT_FRACTION"
            "\tREFERENCE_SQUARED_WEIGHT\tTARGET_AVAILABLE_SQUARED_WEIGHT"
            "\tTARGET_AVAILABLE_SQUARED_WEIGHT_FRACTION\n"
         << std::setprecision(17);
  for (const auto& row : summary.scores) {
    const double row_fraction =
        row.reference_weight_ct
            ? static_cast<double>(row.available_weight_ct) /
                  static_cast<double>(row.reference_weight_ct)
            : 0.0;
    const double l1_fraction =
        row.reference_weight_l1
            ? row.available_weight_l1 / row.reference_weight_l1
            : 0.0;
    const double l2_fraction =
        row.reference_weight_l2_squared
            ? row.available_weight_l2_squared /
                  row.reference_weight_l2_squared
            : 0.0;
    output << row.score_id << '\t' << row.column_name << '\t'
           << row.reference_weight_ct << '\t' << row.available_weight_ct
           << '\t' << row_fraction << '\t' << row.reference_weight_l1 << '\t'
           << row.available_weight_l1 << '\t' << l1_fraction << '\t'
           << row.reference_weight_l2_squared << '\t'
           << row.available_weight_l2_squared << '\t' << l2_fraction << '\n';
  }
  output.close();
  if (!output) throw std::runtime_error("cannot finish " + temporary_output);
  std::filesystem::rename(temporary_output, options.output_path);
  output_cleanup.Release();

  const std::string temporary_metadata = metadata_path + ".tmp";
  RemoveFileOnExit metadata_cleanup(temporary_metadata);
  std::ofstream metadata(temporary_metadata);
  if (!metadata) {
    throw std::runtime_error("cannot write " + temporary_metadata);
  }
  metadata << "{\n"
           << "  \"format\": \"pgensparsescore-fragment-support-v1\",\n"
           << "  \"path\": \""
           << pgensparsescore::JsonEscape(
                  std::filesystem::path(options.output_path).filename().string())
           << "\",\n"
           << "  \"variants\": " << summary.variant_ct << ",\n"
           << "  \"scores\": " << summary.score_ct << ",\n"
           << "  \"reference_weights\": " << summary.reference_weight_ct
           << ",\n"
           << "  \"target_available_weights\": "
           << summary.available_weight_ct << "\n"
           << "}\n";
  metadata.close();
  if (!metadata) {
    throw std::runtime_error("cannot finish " + temporary_metadata);
  }
  std::filesystem::rename(temporary_metadata, metadata_path);
  metadata_cleanup.Release();
  std::cerr << "reported " << summary.available_weight_ct << " of "
            << summary.reference_weight_ct << " reference-projected weights "
            << "across " << summary.score_ct << " scores\n";
  return 0;
}

int FragmentCompileMain(int argc, char** argv) {
  const FragmentCompileOptions options = ParseFragmentCompileOptions(argc, argv);
  if (std::filesystem::exists(options.build.output_path) ||
      std::filesystem::exists(options.build.output_path + ".json")) {
    throw std::runtime_error("score-fragment output already exists: " +
                             options.build.output_path);
  }
  auto progress = MakeProgressReporter(options.progress_jsonl,
                                       options.progress_interval_seconds);
  const auto summary =
      pgensparsescore::CompileScoreFragment(options.build, progress.get());
  std::ofstream metadata(options.build.output_path + ".json");
  if (!metadata) throw std::runtime_error("cannot write fragment metadata");
  metadata << "{\n"
           << "  \"format\": \"pgensparsescore-score-major-fragment-v3\",\n"
           << "  \"path\": \""
           << pgensparsescore::JsonEscape(
                  std::filesystem::path(options.build.output_path)
                      .filename()
                      .string())
           << "\",\n"
           << "  \"variant_index_variants\": "
           << summary.variant_index_variant_ct << ",\n"
           << "  \"tile_size\": " << summary.tile_size << ",\n"
           << "  \"tiles\": " << summary.tile_ct << ",\n"
           << "  \"input_scores\": " << summary.input_score_ct << ",\n"
           << "  \"scores\": " << summary.score_ct << ",\n"
           << "  \"excluded_scores\": " << summary.excluded_score_ct
           << ",\n"
           << "  \"minimum_supported_fraction\": "
           << options.build.minimum_supported_fraction << ",\n"
           << "  \"input_weights\": " << summary.input_weight_ct << ",\n"
           << "  \"catalog_weights\": " << summary.catalog_weight_ct
           << ",\n"
           << "  \"supported_weights\": " << summary.supported_weight_ct
           << ",\n"
           << "  \"weights\": " << summary.weight_ct << ",\n"
           << "  \"zero_weights\": " << summary.zero_weight_ct << ",\n"
           << "  \"excluded_weights\": " << summary.excluded_weight_ct
           << ",\n"
           << "  \"duplicate_weights\": " << summary.duplicate_weight_ct
           << ",\n"
           << "  \"missing_variant_weights\": "
           << summary.missing_variant_weight_ct << ",\n"
           << "  \"missing_frequency_weights\": "
           << summary.missing_frequency_weight_ct << ",\n"
           << "  \"referenced_variants\": "
           << summary.referenced_variant_ct << ",\n"
           << "  \"score_qc_path\": \""
           << pgensparsescore::JsonEscape(
                  std::filesystem::path(options.build.output_path + ".score_qc.tsv")
                      .filename()
                      .string())
           << "\",\n"
           << "  \"variant_bitset_path\": \""
           << pgensparsescore::JsonEscape(
                  std::filesystem::path(options.build.output_path + ".variants.bits")
                      .filename()
                      .string())
           << "\",\n"
           << "  \"variant_bitset_bytes\": "
           << summary.variant_bitset_bytes << ",\n"
           << "  \"bytes\": " << summary.output_bytes << "\n"
           << "}\n";
  metadata.close();
  if (!metadata) throw std::runtime_error("cannot finish fragment metadata");
  std::cerr << "built score fragment with " << summary.score_ct << " of "
            << summary.input_score_ct << " scores and " << summary.weight_ct
            << " retained nonzero weights\n";
  return 0;
}

int VariantBitsMergeMain(int argc, char** argv) {
  const VariantBitsMergeOptions options =
      ParseVariantBitsMergeOptions(argc, argv);
  const auto summary = pgensparsescore::MergeVariantBits(
      options.list_path, options.output_path);
  std::cerr << "merged " << summary.input_ct << " variant bitsets into "
            << summary.referenced_variant_ct << " of " << summary.variant_ct
            << " variants\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "build-variant-index") {
      return VariantIndexMain(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "compile-fragment") {
      return FragmentCompileMain(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "merge-variant-bits") {
      return VariantBitsMergeMain(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "build-support-index") {
      return SupportIndexMain(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "report-fragment-support") {
      return FragmentSupportMain(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "compile") {
      return CompileMain(argc, argv);
    }
    const Options options = ParseOptions(argc, argv);
    const std::filesystem::path output_path(options.out);
    if (!output_path.parent_path().empty()) {
      std::filesystem::create_directories(output_path.parent_path());
    }
    auto progress = MakeProgressReporter(options.progress_jsonl,
                                         options.progress_interval_seconds);
    if (progress) {
      progress->Event("score", "start", {}, {{"output", options.out}});
    }
    const std::vector<pgensparsescore::PfileSpec> inputs =
        options.pfile_list.empty()
            ? std::vector<pgensparsescore::PfileSpec>{
                  {options.pgen, options.pvar, options.psam}}
            : pgensparsescore::ReadPfileList(options.pfile_list);
    if (progress) {
      progress->Event("score", "pfile_inputs_loaded",
                      {{"pgen_inputs", inputs.size()}});
    }
    if (!options.fragment_list.empty()) {
      return RunFragmentScoring(options, inputs, progress.get());
    }
    std::optional<pgensparsescore::CompiledCatalog> compiled_catalog;
    if (!options.compiled_catalog.empty()) {
      compiled_catalog =
          pgensparsescore::ReadCompiledCatalog(options.compiled_catalog,
                                               progress.get());
    }
    std::optional<pgensparsescore::VariantMap> variant_map;
    if (!options.variant_map.empty()) {
      variant_map = pgensparsescore::ReadVariantMap(options.variant_map);
      if (progress) {
        progress->Event("score", "variant_map_loaded",
                        {{"variant_mapping_rows", variant_map->size()}});
      }
    }
    std::unordered_set<std::string> target_variant_ids;
    const bool filter_pvar = compiled_catalog.has_value() ||
                             !options.variant_map.empty();
    if (compiled_catalog) {
      target_variant_ids.reserve(compiled_catalog->variants.size());
      for (const auto& variant : compiled_catalog->variants) {
        if (variant_map) {
          const auto mapping = variant_map->find(variant.source_id);
          if (mapping != variant_map->end()) {
            target_variant_ids.insert(mapping->second);
          }
        } else {
          target_variant_ids.insert(variant.source_id);
        }
      }
    } else if (variant_map) {
      target_variant_ids.reserve(variant_map->size());
      for (const auto& [source_id, target_id] : *variant_map) {
        static_cast<void>(source_id);
        target_variant_ids.insert(target_id);
      }
    }

    auto samples = pgensparsescore::ReadPsam(inputs.front().psam);
    if (progress) {
      progress->Event("score", "samples_loaded",
                      {{"samples", samples.size()}});
    }
    std::vector<pgensparsescore::Variant> all_variants;
    std::vector<uint32_t> input_by_variant;
    std::vector<uint32_t> local_index_by_variant;
    std::vector<uint32_t> pvar_row_counts;
    for (uint32_t input_idx = 0; input_idx < inputs.size(); ++input_idx) {
      if (input_idx) {
        const auto input_samples =
            pgensparsescore::ReadPsam(inputs[input_idx].psam);
        if (!SamplesEqual(samples, input_samples)) {
          throw std::runtime_error(
              "PSAM sample IDs/order differ between PGEN inputs: " +
              inputs[input_idx].psam);
        }
      }
      auto pvar = pgensparsescore::ReadPvar(
          inputs[input_idx].pvar,
          filter_pvar ? &target_variant_ids : nullptr);
      pvar_row_counts.push_back(pvar.row_ct);
      if (all_variants.size() + pvar.variants.size() >
          std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error(
            "combined PVAR inputs exceed the supported variant count");
      }
      for (uint32_t local_idx = 0; local_idx < pvar.variants.size();
           ++local_idx) {
        all_variants.push_back(std::move(pvar.variants[local_idx]));
        input_by_variant.push_back(input_idx);
        local_index_by_variant.push_back(local_idx);
      }
      if (progress) {
        progress->Event(
            "score", "pvar_loaded",
            {{"pgen_input_index", input_idx},
             {"pgen_inputs", inputs.size()},
             {"pvar_rows", pvar.row_ct},
             {"retained_pvar_variants", pvar.variants.size()},
             {"retained_pvar_variants_total", all_variants.size()}},
            {{"pvar", inputs[input_idx].pvar}});
      }
    }
    const uint64_t pvar_variant_ct = all_variants.size();
    auto catalog = compiled_catalog
                       ? pgensparsescore::MaterializeCompiledCatalog(
                             *compiled_catalog, all_variants,
                             variant_map ? &*variant_map : nullptr,
                             progress.get())
                       : pgensparsescore::CompileCatalog(
                             options.manifest, all_variants,
                             variant_map ? &*variant_map : nullptr,
                             progress.get());
    std::optional<pgensparsescore::FrequencyTable> frequencies;
    if (!options.read_freq.empty()) {
      std::unordered_set<std::string> scored_variant_ids;
      scored_variant_ids.reserve(catalog.variants.size());
      for (const auto& variant : catalog.variants) {
        scored_variant_ids.insert(all_variants.at(variant.variant_idx).id);
      }
      frequencies = pgensparsescore::ReadFrequencyTable(
          options.read_freq, &scored_variant_ids);
      if (progress) {
        progress->Event("score", "frequencies_loaded",
                        {{"frequency_rows", frequencies->size()},
                         {"scored_variant_ids", scored_variant_ids.size()}});
      }
    }
    pgensparsescore::ScoreRunStats stats =
        pgensparsescore::ApplyMissingFrequencyPolicy(
            &catalog, all_variants,
            frequencies ? &*frequencies : nullptr,
            options.missing_frequency_policy);
    auto catalog_by_input =
        PartitionCatalog(&catalog, input_by_variant, local_index_by_variant,
                         inputs.size());
    std::vector<std::vector<pgensparsescore::Variant>> variants_by_input(
        inputs.size());
    for (size_t variant_idx = 0; variant_idx < all_variants.size();
         ++variant_idx) {
      variants_by_input[input_by_variant[variant_idx]].push_back(
          std::move(all_variants[variant_idx]));
    }
    all_variants.clear();
    all_variants.shrink_to_fit();
    input_by_variant.clear();
    input_by_variant.shrink_to_fit();
    local_index_by_variant.clear();
    local_index_by_variant.shrink_to_fit();
    const auto output_format = ParseScoreOutputFormat(options.output_format);
    const std::string working_path =
        output_format == ScoreOutputFormat::kScoreMajorBinary
            ? ScoreOutputPath(options.out, output_format) + ".tmp"
            : options.out + ".work.score-major.bin";
    RemoveFileOnExit remove_working(working_path);
    uint64_t working_matrix_byte_ct = 0;
    std::vector<std::string> pgen_storage_modes(inputs.size());
    {
      pgensparsescore::MappedMatrix matrix(
          working_path, catalog.scores.size(), samples.size());
      working_matrix_byte_ct = matrix.byte_ct();
      if (progress) {
        progress->Event("score", "working_matrix_ready",
                        {{"working_matrix_bytes", working_matrix_byte_ct},
                         {"sample_rows", samples.size()},
                         {"score_columns", catalog.scores.size()}});
      }

      auto process = [&](size_t input_idx,
                         const pgensparsescore::PfileSpec& input,
                         const std::vector<pgensparsescore::Variant>& variants,
                         const pgensparsescore::Catalog& input_catalog,
                         uint32_t pvar_row_ct) {
        pgensparsescore::PgenDosageReader reader(input.pgen);
        pgen_storage_modes[input_idx] = reader.storage_mode_name();
        if (reader.variant_ct() != pvar_row_ct) {
          throw std::runtime_error("PGEN/PVAR variant-count mismatch: " +
                                   input.pgen);
        }
        if (reader.sample_ct() != samples.size()) {
          throw std::runtime_error("PGEN/PSAM sample-count mismatch: " +
                                   input.pgen);
        }
        if (progress) {
          progress->Event(
              "score", "pgen_start",
              {{"pgen_input_index", input_idx},
               {"pgen_inputs", inputs.size()},
               {"variants_to_score", input_catalog.variants.size()}},
              {{"pgen", input.pgen},
               {"storage_mode", reader.storage_mode_name()}});
        }
        const auto input_stats = pgensparsescore::ScoreCatalog(
            input_catalog, variants, frequencies ? &*frequencies : nullptr,
            options.missing_frequency_policy, &reader, &matrix,
            progress.get());
        AddStats(input_stats, &stats);
        if (progress) {
          progress->Event(
              "score", "pgen_scored",
              {{"pgen_input_index", input_idx},
               {"pgen_inputs", inputs.size()},
               {"scored_variants_total", stats.variant_ct},
               {"weight_edges_total", stats.edge_ct},
               {"sparse_variants_total", stats.sparse_variant_ct},
               {"dense_variants_total", stats.dense_variant_ct},
               {"sparse_weight_edges_total", stats.sparse_edge_ct},
               {"dense_weight_edges_total", stats.dense_edge_ct},
               {"sparse_score_updates_total", stats.sparse_update_ct},
               {"dense_score_updates_total", stats.dense_update_ct}});
        }
      };

      for (size_t input_idx = 0; input_idx < inputs.size(); ++input_idx) {
        process(input_idx, inputs[input_idx], variants_by_input[input_idx],
                catalog_by_input[input_idx], pvar_row_counts[input_idx]);
      }
      if (output_format == ScoreOutputFormat::kWideTsv) {
        WriteWideScores(options.out, samples, catalog, matrix, progress.get());
      } else {
        FinishBinaryScores(options.out, samples, catalog, &matrix,
                           working_path, &remove_working, progress.get());
      }
    }
    if (output_format == ScoreOutputFormat::kWideTsv) {
      remove_working.RemoveNow();
    }
    WriteMetadata(options.out, samples.size(), samples.front().fid.has_value(),
                  inputs.size(), frequencies ? frequencies->size() : 0,
                  variant_map ? variant_map->size() : 0,
                  pvar_variant_ct, pgen_storage_modes,
                  compiled_catalog ? "compiled" : "manifest",
                  options.missing_frequency_policy, working_matrix_byte_ct,
                  catalog, stats, 0, 0, 1,
                  pgensparsescore::DenseScoringKernel::kDirect, 0,
                  output_format);
    if (progress) {
      progress->Event(
          "score", "complete",
          {{"sample_rows", samples.size()},
           {"score_columns", catalog.scores.size()},
           {"pgen_inputs", inputs.size()},
           {"scored_variants", stats.variant_ct},
           {"weight_edges", stats.edge_ct},
           {"sparse_variants", stats.sparse_variant_ct},
           {"dense_variants", stats.dense_variant_ct},
           {"sparse_weight_edges", stats.sparse_edge_ct},
           {"dense_weight_edges", stats.dense_edge_ct},
           {"sparse_score_updates", stats.sparse_update_ct},
           {"dense_score_updates", stats.dense_update_ct},
           {"imputed_values", stats.imputed_value_ct}},
          {{"output_format", ScoreOutputFormatName(output_format)}});
    }
    std::cerr << "wrote " << catalog.scores.size()
              << " named score columns for " << samples.size()
              << " sample rows from " << inputs.size() << " PGEN input(s); "
              << stats.sparse_variant_ct << " sparse and "
              << stats.dense_variant_ct << " dense variant decodes\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "pgensparsescore: " << error.what() << '\n';
    PrintUsage(std::cerr);
    return 1;
  }
}
