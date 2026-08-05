// SPDX-License-Identifier: GPL-3.0-only
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unistd.h>
#include <vector>
#include <zstd.h>

#include "catalog.h"
#include "compiled_catalog.h"
#include "frequency.h"
#include "mapped_matrix.h"
#include "pfile.h"
#include "scorer.h"

namespace {

void ExpectNear(double actual, double expected, const std::string& message) {
  if (std::abs(actual - expected) > 1e-12) {
    throw std::runtime_error(message + ": got " + std::to_string(actual) +
                             ", expected " + std::to_string(expected));
  }
}

std::filesystem::path TempPath(const std::string& suffix) {
  return std::filesystem::temp_directory_path() /
         ("pgensparsescore-test-" + std::to_string(::getpid()) + suffix);
}

void TestDenseKernel() {
  const auto path = TempPath("-dense.bin");
  {
    pgensparsescore::MappedMatrix matrix(path.string(), 2, 4);
    const double dosages[] = {0.0, 1.0, 2.0, 0.5};
    const std::vector<pgensparsescore::Edge> edges{{0, 2.0}, {1, -1.0}};
    pgensparsescore::ApplyDenseDosage(dosages, 4, edges, &matrix);
    const double expected0[] = {0.0, 2.0, 4.0, 1.0};
    const double expected1[] = {0.0, -1.0, -2.0, -0.5};
    for (uint32_t idx = 0; idx < 4; ++idx) {
      ExpectNear(matrix.Row(0)[idx], expected0[idx], "dense row 0");
      ExpectNear(matrix.Row(1)[idx], expected1[idx], "dense row 1");
    }
  }
  std::filesystem::remove(path);
}

void TestSparseKernel() {
  const auto path = TempPath("-sparse.bin");
  {
    pgensparsescore::MappedMatrix matrix(path.string(), 1, 5);
    const uint32_t sample_ids[] = {1, 3, 4};
    const uint16_t values[] = {16384, UINT16_MAX, 32768};
    const std::vector<pgensparsescore::Edge> edges{{0, 3.0}};
    std::vector<double> baselines{0.0};
    pgensparsescore::ApplySparseDosage(0.0, 1.0, sample_ids, values, 3,
                                       edges, &baselines, &matrix);
    ExpectNear(baselines[0], 0.0, "sparse baseline");
    const double expected[] = {0.0, 3.0, 0.0, 3.0, 6.0};
    for (uint32_t idx = 0; idx < 5; ++idx) {
      ExpectNear(matrix.Row(0)[idx], expected[idx], "sparse row");
    }
  }
  std::filesystem::remove(path);
}

void TestCatalogOrientation() {
  const auto directory = TempPath("-catalog");
  std::filesystem::create_directories(directory);
  const auto weight_path = directory / "weights.tsv";
  const auto manifest_path = directory / "manifest.tsv";
  {
    std::ofstream output(weight_path);
    output << "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"
           << "v1\tG\tA\t2\n"
           << "v2\tC\tT\t3\n"
           << "absent\tA\tG\t4\n";
  }
  {
    std::ofstream output(manifest_path);
    output << "SCORE\tPATH\nscore1\tweights.tsv\n";
  }
  const std::vector<pgensparsescore::Variant> variants{
      {"1", "v1", "A", "G"}, {"1", "v2", "C", "T"}};
  const auto catalog =
      pgensparsescore::CompileCatalog(manifest_path.string(), variants);
  if (catalog.scores.size() != 1 || catalog.variants.size() != 2) {
    throw std::runtime_error("catalog dimensions are wrong");
  }
  ExpectNear(catalog.intercepts[0], 6.0, "REF effect intercept");
  ExpectNear(catalog.variants[0].edges[0].beta_alt, 2.0, "ALT beta");
  ExpectNear(catalog.variants[1].edges[0].beta_alt, -3.0, "REF beta");
  if (catalog.scores[0].missing_variant_ct != 1) {
    throw std::runtime_error("missing variant count is wrong");
  }
  std::filesystem::remove_all(directory);
}

void TestRepeatedWeightRowsAreRetained() {
  const auto directory = TempPath("-duplicate-weight");
  std::filesystem::create_directories(directory);
  const auto weight_path = directory / "weights.tsv";
  const auto manifest_path = directory / "manifest.tsv";
  const auto catalog_path = directory / "catalog.bin";
  {
    std::ofstream output(weight_path);
    output << "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"
           << "v1\tG\tA\t2\n"
           << "v1\tG\tA\t3\n";
  }
  {
    std::ofstream output(manifest_path);
    output << "SCORE\tPATH\nscore1\tweights.tsv\n";
  }
  const auto compiled =
      pgensparsescore::CompileSourceCatalog(manifest_path.string());
  if (compiled.weight_ct != 2 || compiled.variants.size() != 1 ||
      compiled.scores[0].input_weight_ct != 2 ||
      compiled.scores[0].duplicate_weight_ct != 1 ||
      compiled.scores[0].catalog_weight_ct != 2 ||
      compiled.variants[0].weights.size() != 2) {
    throw std::runtime_error("repeated weight rows were not retained");
  }
  pgensparsescore::WriteCompiledCatalog(catalog_path.string(), compiled);
  const auto round_trip =
      pgensparsescore::ReadCompiledCatalog(catalog_path.string());
  if (round_trip.scores[0].duplicate_weight_ct != 1 ||
      round_trip.weight_ct != 2 || round_trip.variants[0].weights.size() != 2) {
    throw std::runtime_error("duplicate count did not survive round trip");
  }
  const std::vector<pgensparsescore::Variant> variants{
      {"1", "v1", "A", "G"}};
  const auto materialized =
      pgensparsescore::MaterializeCompiledCatalog(round_trip, variants, nullptr);
  if (materialized.scores[0].matched_weight_ct != 2 ||
      materialized.variants[0].edges.size() != 2) {
    throw std::runtime_error("repeated rows were not materialized");
  }
  ExpectNear(materialized.variants[0].edges[0].beta_alt +
                 materialized.variants[0].edges[1].beta_alt,
             5.0, "repeated row coefficient sum");
  std::filesystem::remove_all(directory);
}

void TestVariantMapping() {
  const auto directory = TempPath("-variant-map");
  std::filesystem::create_directories(directory);
  const auto weight_path = directory / "weights.tsv";
  const auto manifest_path = directory / "manifest.tsv";
  const auto mapping_path = directory / "mapping.tsv";
  {
    std::ofstream output(weight_path);
    output << "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"
           << "source-v1\tG\tA\t2\n"
           << "unmapped\tA\tG\t4\n";
  }
  {
    std::ofstream output(manifest_path);
    output << "SCORE\tPATH\nscore1\tweights.tsv\n";
  }
  {
    std::ofstream output(mapping_path);
    output << "SOURCE_ID\tTARGET_ID\nsource-v1\ttarget-v1\n";
  }
  const auto mapping =
      pgensparsescore::ReadVariantMap(mapping_path.string());
  const std::vector<pgensparsescore::Variant> variants{
      {"1", "target-v1", "A", "G"}};
  const auto catalog = pgensparsescore::CompileCatalog(
      manifest_path.string(), variants, &mapping);
  if (catalog.variants.size() != 1 ||
      catalog.scores[0].matched_weight_ct != 1 ||
      catalog.scores[0].missing_variant_ct != 1) {
    throw std::runtime_error("variant-map matching/QC counts are wrong");
  }
  ExpectNear(catalog.variants[0].edges[0].beta_alt, 2.0,
             "mapped ALT beta");
  std::filesystem::remove_all(directory);
}

void TestFilteredPvarKeepsPgenIndexes() {
  const auto directory = TempPath("-filtered-pvar");
  std::filesystem::create_directories(directory);
  const auto pvar_path = directory / "source.pvar";
  {
    std::ofstream output(pvar_path);
    output << "#CHROM\tPOS\tID\tREF\tALT\n"
           << "1\t100\tv1\tA\tG\n"
           << "1\t200\tv2\tC\tT\n"
           << "2\t300\tv3\tG\tA\n";
  }
  const std::unordered_set<std::string> included{"v2"};
  const auto pvar =
      pgensparsescore::ReadPvar(pvar_path.string(), &included);
  if (pvar.row_ct != 3 || pvar.variants.size() != 1 ||
      pvar.variants[0].id != "v2" ||
      pvar.variants[0].pgen_variant_idx != 1) {
    throw std::runtime_error("filtered PVAR did not retain source row index");
  }
  std::filesystem::remove_all(directory);
}

void WriteZstd(const std::filesystem::path& path, const std::string& contents) {
  std::vector<char> compressed(ZSTD_compressBound(contents.size()));
  const size_t compressed_size =
      ZSTD_compress(compressed.data(), compressed.size(), contents.data(),
                    contents.size(), 1);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error("cannot create zstd test input");
  }
  std::ofstream output(path, std::ios::binary);
  output.write(compressed.data(), static_cast<std::streamsize>(compressed_size));
}

void TestFrequencyParsing() {
  const auto directory = TempPath("-frequency");
  std::filesystem::create_directories(directory);
  const std::string contents =
      "#CHROM\tID\tREF\tALT\tREF_CT\tALT_CTS\tOBS_CT\n"
      "1\tv1\tA\tG\t5\t1\t6\n"
      "2\tv2\tC\tT\t0\tT=8\t8\n"
      "3\tall-missing\tG\tA\t0\t0\t0\n";
  const auto plain_path = directory / "frequency.acount";
  const auto zstd_path = directory / "frequency.acount.zst";
  {
    std::ofstream output(plain_path);
    output << contents;
  }
  WriteZstd(zstd_path, contents);

  for (const auto& path : {plain_path, zstd_path}) {
    const auto frequencies =
        pgensparsescore::ReadFrequencyTable(path.string());
    if (frequencies.size() != 2 || frequencies.at("v1").ref != "A" ||
        frequencies.at("v1").alt != "G") {
      throw std::runtime_error("frequency metadata is wrong");
    }
    ExpectNear(frequencies.at("v1").alt_dosage_mean, 1.0 / 3.0,
               "ALT count dosage mean");
    ExpectNear(frequencies.at("v2").alt_dosage_mean, 2.0,
               "allele=value ALT count dosage mean");
  }
  std::filesystem::remove_all(directory);
}

void TestPfileList() {
  const auto directory = TempPath("-pfile-list");
  std::filesystem::create_directories(directory / "inputs");
  const auto list_path = directory / "pfiles.tsv";
  {
    std::ofstream output(list_path);
    output << "PGEN\tPVAR\tPSAM\n"
           << "inputs/chr1.pgen\tinputs/chr1.pvar.zst\tinputs/chr1.psam\n";
  }
  const auto inputs = pgensparsescore::ReadPfileList(list_path.string());
  if (inputs.size() != 1 ||
      inputs[0].pgen != (directory / "inputs/chr1.pgen").string() ||
      inputs[0].pvar != (directory / "inputs/chr1.pvar.zst").string() ||
      inputs[0].psam != (directory / "inputs/chr1.psam").string()) {
    throw std::runtime_error("PGEN list path resolution is wrong");
  }
  std::filesystem::remove_all(directory);
}

}  // namespace

int main() {
  try {
    TestDenseKernel();
    TestSparseKernel();
    TestCatalogOrientation();
    TestRepeatedWeightRowsAreRetained();
    TestVariantMapping();
    TestFilteredPvarKeepsPgenIndexes();
    TestFrequencyParsing();
    TestPfileList();
    std::cout << "all unit tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
