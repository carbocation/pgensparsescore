// SPDX-License-Identifier: GPL-3.0-only
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "catalog.h"
#include "mapped_matrix.h"
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

}  // namespace

int main() {
  try {
    TestDenseKernel();
    TestSparseKernel();
    TestCatalogOrientation();
    std::cout << "all unit tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
