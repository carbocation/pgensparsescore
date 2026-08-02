// SPDX-License-Identifier: GPL-3.0-only
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

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

void WriteMetadata(const std::string& prefix,
                   const std::vector<pgensparsescore::Sample>& samples,
                   const pgensparsescore::Catalog& catalog,
                   const pgensparsescore::ScoreRunStats& stats,
                   uint64_t matrix_byte_ct) {
  {
    std::ofstream output(prefix + ".samples.tsv");
    if (!output) throw std::runtime_error("cannot write sample metadata");
    output << "INDEX\tFID\tIID\n";
    for (uint32_t idx = 0; idx < samples.size(); ++idx) {
      output << idx << '\t' << samples[idx].fid << '\t' << samples[idx].iid
             << '\n';
    }
  }
  {
    std::ofstream output(prefix + ".scores.tsv");
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
    const std::filesystem::path binary(prefix + ".scores.bin");
    output << "{\n"
           << "  \"format\": \"pgensparsescore-matrix-v1\",\n"
           << "  \"path\": \""
           << pgensparsescore::JsonEscape(binary.filename().string()) << "\",\n"
           << "  \"dtype\": \"<f8\",\n"
           << "  \"order\": \"C\",\n"
           << "  \"shape\": [" << catalog.scores.size() << ", "
           << samples.size() << "],\n"
           << "  \"matrix_bytes\": " << matrix_byte_ct << ",\n"
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
    pgensparsescore::MappedMatrix matrix(options.out + ".scores.bin",
                                         catalog.scores.size(), samples.size());
    const auto stats =
        pgensparsescore::ScoreCatalog(catalog, &reader, &matrix);
    WriteMetadata(options.out, samples, catalog, stats, matrix.byte_ct());
    std::cerr << "wrote " << catalog.scores.size() << " scores for "
              << samples.size() << " samples; " << stats.sparse_variant_ct
              << " sparse and " << stats.dense_variant_ct
              << " dense variant decodes\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "pgensparsescore: " << error.what() << '\n';
    PrintUsage(std::cerr);
    return 1;
  }
}
