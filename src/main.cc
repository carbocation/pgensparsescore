// SPDX-License-Identifier: GPL-3.0-only
#include <algorithm>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "catalog.h"
#include "io.h"
#include "mapped_matrix.h"
#include "pgen_reader.h"
#include "scorer.h"

namespace {

struct Options {
  std::string pgen;
  std::string pvar;
  std::string psam;
  std::string manifest;
  std::string out;
};

void PrintUsage(std::ostream& stream) {
  stream <<
      "Usage: pgensparsescore --pgen FILE --pvar FILE --psam FILE\n"
      "                       --manifest FILE --out PREFIX\n";
}

Options ParseOptions(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--help") {
    PrintUsage(std::cout);
    std::exit(0);
  }
  std::unordered_map<std::string, std::string*> destinations;
  Options options;
  destinations["--pgen"] = &options.pgen;
  destinations["--pvar"] = &options.pvar;
  destinations["--psam"] = &options.psam;
  destinations["--manifest"] = &options.manifest;
  destinations["--out"] = &options.out;
  for (int idx = 1; idx < argc; ++idx) {
    const std::string key(argv[idx]);
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
  for (const auto& [name, destination] : destinations) {
    if (destination->empty()) {
      throw std::runtime_error("required argument is missing: " + name);
    }
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
  pgensparsescore::GzipWriter output(temporary_path);

  std::string line = "FID\tIID";
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
      line += sample.fid;
      line.push_back('\t');
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
}

void WriteMetadata(const std::string& prefix, uint32_t sample_ct,
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
           << "  \"sample_rows\": " << sample_ct << ",\n"
           << "  \"score_columns\": " << catalog.scores.size() << ",\n"
           << "  \"working_matrix_bytes\": " << working_matrix_byte_ct
           << ",\n"
           << "  \"scored_variants\": " << stats.variant_ct << ",\n"
           << "  \"weight_edges\": " << stats.edge_ct << ",\n"
           << "  \"sparse_variants\": " << stats.sparse_variant_ct << ",\n"
           << "  \"dense_variants\": " << stats.dense_variant_ct << ",\n"
           << "  \"imputed_values\": " << stats.imputed_value_ct << "\n"
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
    const auto variants = pgensparsescore::ReadPvar(options.pvar);
    const auto samples = pgensparsescore::ReadPsam(options.psam);
    pgensparsescore::PgenDosageReader reader(options.pgen);
    if (reader.variant_ct() != variants.size()) {
      throw std::runtime_error("PGEN/PVAR variant-count mismatch");
    }
    if (reader.sample_ct() != samples.size()) {
      throw std::runtime_error("PGEN/PSAM sample-count mismatch");
    }
    const auto catalog =
        pgensparsescore::CompileCatalog(options.manifest, variants);

    const std::string working_path = options.out + ".work.score-major.bin";
    pgensparsescore::ScoreRunStats stats;
    uint64_t working_matrix_byte_ct = 0;
    {
      pgensparsescore::MappedMatrix matrix(
          working_path, catalog.scores.size(), samples.size());
      working_matrix_byte_ct = matrix.byte_ct();
      stats = pgensparsescore::ScoreCatalog(catalog, &reader, &matrix);
      WriteWideScores(options.out, samples, catalog, matrix);
    }
    if (!std::filesystem::remove(working_path)) {
      throw std::runtime_error("cannot remove working score matrix " +
                               working_path);
    }
    WriteMetadata(options.out, samples.size(), working_matrix_byte_ct, catalog,
                  stats);
    std::cerr << "wrote " << catalog.scores.size() << " named score columns for "
              << samples.size() << " sample rows; " << stats.sparse_variant_ct
              << " sparse and " << stats.dense_variant_ct
              << " dense variant decodes\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "pgensparsescore: " << error.what() << '\n';
    PrintUsage(std::cerr);
    return 1;
  }
}
