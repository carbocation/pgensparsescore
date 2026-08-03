// SPDX-License-Identifier: GPL-3.0-only
#include <algorithm>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "catalog.h"
#include "frequency.h"
#include "io.h"
#include "mapped_matrix.h"
#include "pfile.h"
#include "pgen_reader.h"
#include "scorer.h"

namespace {

struct Options {
  std::string pgen;
  std::string pvar;
  std::string psam;
  std::string pfile_list;
  std::string manifest;
  std::string variant_map;
  std::string read_freq;
  std::string out;
  bool error_on_missing_freq = false;
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
      "  pgensparsescore --pgen FILE --pvar FILE --psam FILE \\\n"
      "                     --manifest FILE [--variant-map FILE] \\\n"
      "                     [--read-freq FILE] \\\n"
      "                     [--error-on-missing-freq] --out PREFIX\n"
      "  pgensparsescore --pfile-list FILE --manifest FILE \\\n"
      "                     [--variant-map FILE] \\\n"
      "                     [--read-freq FILE] \\\n"
      "                     [--error-on-missing-freq] --out PREFIX\n";
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
      {"--variant-map", &options.variant_map},
      {"--read-freq", &options.read_freq},
      {"--out", &options.out},
  };
  for (int idx = 1; idx < argc; ++idx) {
    const std::string key(argv[idx]);
    if (key == "--error-on-missing-freq") {
      if (options.error_on_missing_freq) {
        throw std::runtime_error("argument supplied twice: " + key);
      }
      options.error_on_missing_freq = true;
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
  if (options.manifest.empty() || options.out.empty()) {
    throw std::runtime_error("--manifest and --out are required");
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
  if (options.error_on_missing_freq && options.read_freq.empty()) {
    throw std::runtime_error(
        "--error-on-missing-freq requires --read-freq");
  }
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
                     const pgensparsescore::MappedMatrix& matrix) {
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
  }
  output.Close();
  std::filesystem::rename(temporary_path, final_path);
  remove_temporary.Release();
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
  output->imputed_value_ct += input.imputed_value_ct;
  output->external_frequency_variant_ct += input.external_frequency_variant_ct;
  output->cohort_frequency_variant_ct += input.cohort_frequency_variant_ct;
  output->missing_frequency_variant_ct += input.missing_frequency_variant_ct;
}

void WriteMetadata(const std::string& prefix, uint32_t sample_ct, bool has_fid,
                   uint32_t pgen_ct, uint64_t frequency_row_ct,
                   uint64_t variant_mapping_row_ct,
                   bool error_on_missing_frequency,
                   uint64_t working_matrix_byte_ct,
                   const pgensparsescore::Catalog& catalog,
                   const pgensparsescore::ScoreRunStats& stats) {
  {
    std::ofstream output(prefix + ".score-metadata.tsv");
    if (!output) throw std::runtime_error("cannot write score metadata");
    output << "INDEX\tSCORE\tINPUT_WEIGHTS\tMATCHED_WEIGHTS\tMISSING_VARIANTS"
              "\tALT_EFFECTS\tREF_EFFECTS\tREF_INTERCEPT\n";
    for (uint32_t idx = 0; idx < catalog.scores.size(); ++idx) {
      const auto& score = catalog.scores[idx];
      output << idx << '\t' << score.id << '\t' << score.input_weight_ct << '\t'
             << score.matched_weight_ct << '\t' << score.missing_variant_ct
             << '\t' << score.alt_effect_ct << '\t' << score.ref_effect_ct
             << '\t' << score.ref_effect_intercept << '\n';
    }
  }
  {
    std::ofstream output(prefix + ".json");
    if (!output) throw std::runtime_error("cannot write JSON metadata");
    const std::filesystem::path scores(prefix + ".scores.tsv.gz");
    output << "{\n"
           << "  \"format\": \"pgensparsescore-wide-tsv-v1\",\n"
           << "  \"path\": \""
           << pgensparsescore::JsonEscape(scores.filename().string()) << "\",\n"
           << "  \"sample_id_columns\": "
           << (has_fid ? "[\"FID\", \"IID\"]" : "[\"IID\"]") << ",\n"
           << "  \"sample_rows\": " << sample_ct << ",\n"
           << "  \"score_columns\": " << catalog.scores.size() << ",\n"
           << "  \"pgen_inputs\": " << pgen_ct << ",\n"
           << "  \"frequency_rows\": " << frequency_row_ct << ",\n"
           << "  \"variant_mapping_rows\": " << variant_mapping_row_ct
           << ",\n"
           << "  \"error_on_missing_frequency\": "
           << (error_on_missing_frequency ? "true" : "false") << ",\n"
           << "  \"working_matrix_bytes\": " << working_matrix_byte_ct
           << ",\n"
           << "  \"scored_variants\": " << stats.variant_ct << ",\n"
           << "  \"weight_edges\": " << stats.edge_ct << ",\n"
           << "  \"sparse_variants\": " << stats.sparse_variant_ct << ",\n"
           << "  \"dense_variants\": " << stats.dense_variant_ct << ",\n"
           << "  \"imputed_values\": " << stats.imputed_value_ct << ",\n"
           << "  \"external_frequency_variants\": "
           << stats.external_frequency_variant_ct << ",\n"
           << "  \"cohort_frequency_variants\": "
           << stats.cohort_frequency_variant_ct << ",\n"
           << "  \"missing_frequency_variants\": "
           << stats.missing_frequency_variant_ct << "\n"
           << "}\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = ParseOptions(argc, argv);
    const std::filesystem::path output_path(options.out);
    if (!output_path.parent_path().empty()) {
      std::filesystem::create_directories(output_path.parent_path());
    }
    const std::vector<pgensparsescore::PfileSpec> inputs =
        options.pfile_list.empty()
            ? std::vector<pgensparsescore::PfileSpec>{
                  {options.pgen, options.pvar, options.psam}}
            : pgensparsescore::ReadPfileList(options.pfile_list);
    std::optional<pgensparsescore::FrequencyTable> frequencies;
    if (!options.read_freq.empty()) {
      frequencies = pgensparsescore::ReadFrequencyTable(options.read_freq);
    }
    std::optional<pgensparsescore::VariantMap> variant_map;
    if (!options.variant_map.empty()) {
      variant_map = pgensparsescore::ReadVariantMap(options.variant_map);
    }

    auto samples = pgensparsescore::ReadPsam(inputs.front().psam);
    std::vector<pgensparsescore::Variant> all_variants;
    std::vector<uint32_t> input_by_variant;
    std::vector<uint32_t> local_index_by_variant;
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
      auto variants = pgensparsescore::ReadPvar(inputs[input_idx].pvar);
      if (all_variants.size() + variants.size() >
          std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error(
            "combined PVAR inputs exceed the supported variant count");
      }
      for (uint32_t local_idx = 0; local_idx < variants.size(); ++local_idx) {
        all_variants.push_back(std::move(variants[local_idx]));
        input_by_variant.push_back(input_idx);
        local_index_by_variant.push_back(local_idx);
      }
    }
    auto catalog = pgensparsescore::CompileCatalog(
        options.manifest, all_variants,
        variant_map ? &*variant_map : nullptr);
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
    pgensparsescore::ScoreRunStats stats;

    const std::string working_path = options.out + ".work.score-major.bin";
    RemoveFileOnExit remove_working(working_path);
    uint64_t working_matrix_byte_ct = 0;
    {
      pgensparsescore::MappedMatrix matrix(
          working_path, catalog.scores.size(), samples.size());
      working_matrix_byte_ct = matrix.byte_ct();

      auto process = [&](const pgensparsescore::PfileSpec& input,
                         const std::vector<pgensparsescore::Variant>& variants,
                         const pgensparsescore::Catalog& input_catalog) {
        pgensparsescore::PgenDosageReader reader(input.pgen);
        if (reader.variant_ct() != variants.size()) {
          throw std::runtime_error("PGEN/PVAR variant-count mismatch: " +
                                   input.pgen);
        }
        if (reader.sample_ct() != samples.size()) {
          throw std::runtime_error("PGEN/PSAM sample-count mismatch: " +
                                   input.pgen);
        }
        const auto input_stats = pgensparsescore::ScoreCatalog(
            input_catalog, variants, frequencies ? &*frequencies : nullptr,
            options.error_on_missing_freq, &reader, &matrix);
        AddStats(input_stats, &stats);
      };

      for (size_t input_idx = 0; input_idx < inputs.size(); ++input_idx) {
        process(inputs[input_idx], variants_by_input[input_idx],
                catalog_by_input[input_idx]);
      }
      WriteWideScores(options.out, samples, catalog, matrix);
    }
    remove_working.RemoveNow();
    WriteMetadata(options.out, samples.size(), samples.front().fid.has_value(),
                  inputs.size(), frequencies ? frequencies->size() : 0,
                  variant_map ? variant_map->size() : 0,
                  options.error_on_missing_freq, working_matrix_byte_ct,
                  catalog, stats);
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
