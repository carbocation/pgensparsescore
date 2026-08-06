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
#include "fragment_scorer.h"
#include "mapped_matrix.h"
#include "pgen_rans_container.h"
#include "pgen_rans_hybrid.h"
#include "pfile.h"
#include "pgen_reader.h"
#include "scorer.h"
#include "score_fragment.h"
#include "support_index.h"
#include "variant_index.h"

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

void TestConditionalRansPgenReader() {
  constexpr uint32_t sample_ct = 33;
  constexpr uint32_t variant_ct = 2;
  const auto path = TempPath("-conditional-rans.pgen");
  std::vector<uint64_t> sparse(pgen_rans::PackedWordCt(sample_ct), 0);
  pgen_rans::SetPackedGenotype(sparse.data(), 3, 1);
  pgen_rans::SetPackedGenotype(sparse.data(), 17, 3);
  std::vector<uint64_t> dense(pgen_rans::PackedWordCt(sample_ct), 0);
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    pgen_rans::SetPackedGenotype(
        dense.data(), sample_idx,
        static_cast<uint8_t>((sample_idx + 1) % 3));
  }
  pgen_rans::SetPackedGenotype(dense.data(), 16, 3);

  pgen_rans::EncodedBlock block;
  block.first_variant = 0;
  block.records.resize(variant_ct);
  std::string error;
  if (!pgen_rans::EncodeSparsePredictorRecord(
          sparse.data(), nullptr, nullptr, sample_ct,
          pgen_rans::RecordMode::kMarginal, 0, 0, &block.records[0],
          &error) ||
      !pgen_rans::EncodeRawRecord(
          dense.data(), sample_ct, &block.records[1], &error)) {
    throw std::runtime_error("cannot encode conditional-rANS fixture: " +
                             error);
  }
  pgen_rans::ContainerWriter writer;
  const pgen_rans::ContainerParams params(
      sample_ct, variant_ct, variant_ct, 1, 4, 12, 1, 1, 2);
  if (!writer.Open(path.string(), params, {}, &error) ||
      !writer.WriteBlock(block, &error) || !writer.Close(&error)) {
    throw std::runtime_error("cannot write conditional-rANS fixture: " +
                             error);
  }

  {
    pgensparsescore::PgenDosageReader reader(path.string());
    if (reader.sample_ct() != sample_ct || reader.variant_ct() != variant_ct ||
        std::string(reader.storage_mode_name()) != "conditional-rans") {
      throw std::runtime_error("conditional-rANS PGEN metadata is wrong");
    }
    const auto sparse_view = reader.Read(0, 0.25);
    if (!sparse_view.sparse || sparse_view.common != 0.0 ||
        sparse_view.sparse_value_ct != 2 || sparse_view.missing_ct != 1 ||
        sparse_view.sparse_sample_ids[0] != 3 ||
        sparse_view.sparse_sample_ids[1] != 17 ||
        sparse_view.sparse_dosage16[0] != 16384 ||
        sparse_view.sparse_dosage16[1] != UINT16_MAX) {
      throw std::runtime_error("conditional-rANS sparse dosage is wrong");
    }
    const auto dense_view = reader.Read(1, 0.75);
    if (dense_view.sparse || dense_view.missing_ct != 1) {
      throw std::runtime_error("conditional-rANS dense dosage is wrong");
    }
    ExpectNear(dense_view.dense_values[16], 0.75,
               "conditional-rANS missing-value imputation");
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
    output << "SCORE_ID\tCOLUMN_NAME\tDISPLAY_NAME\tPATH\n"
           << "source:score1\tsource__score1\tsource__score1__Trait"
              "\tweights.tsv\n";
  }
  const std::vector<pgensparsescore::Variant> variants{
      {"1", "v1", "A", "G"}, {"1", "v2", "C", "T"}};
  const auto catalog =
      pgensparsescore::CompileCatalog(manifest_path.string(), variants);
  if (catalog.scores.size() != 1 || catalog.variants.size() != 2 ||
      catalog.scores[0].id != "source__score1") {
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

void TestVariantIndex() {
  const auto directory = TempPath("-variant-index");
  std::filesystem::create_directories(directory);
  const auto variants_path = directory / "variants.tsv";
  const auto index_path = directory / "variants.index.bin";
  {
    std::ofstream output(variants_path);
    output << "AOU_ID\tCANONICAL_KEY\tREF\tALT\n"
           << "DRAGEN:chr1:100:A:G\tchr1:100:A:G\tA\tG\n"
           << "same-id\tsame-id\tC\tT\n"
           << "DRAGEN:chr2:300:AT:A\tchr2:300:AT:A\tAT\tA\n";
  }
  pgensparsescore::VariantIndexBuildOptions options;
  options.input_path = variants_path.string();
  options.source_id_column = "AOU_ID";
  options.target_id_column = "CANONICAL_KEY";
  options.block_size = 2;
  options.output_path = index_path.string();
  pgensparsescore::BuildVariantIndex(options);
  {
    pgensparsescore::VariantIndex index(index_path.string());
    if (index.variant_ct() != 3 || index.alias_ct() != 5 ||
        index.block_size() != 2 || index.block_ct() != 2 ||
        index.Lookup("DRAGEN:chr1:100:A:G") != 0 ||
        index.Lookup("chr1:100:A:G") != 0 ||
        index.Lookup("same-id") != 1 ||
        index.Lookup("chr2:300:AT:A") != 2 ||
        index.Lookup("absent").has_value() || index.ref(2) != "AT" ||
        index.alt(2) != "A") {
      throw std::runtime_error("variant index lookup or metadata is wrong");
    }
  }

  const auto prefixed_path = directory / "prefixed.tsv";
  const auto prefixed_index_path = directory / "prefixed.index.bin";
  {
    std::ofstream output(prefixed_path);
    output << "SOURCE_ID\tTARGET_ID\tREF\tALT\n"
           << "DRAGEN:chr1:400:G:C\tchr1:400:G:C\tG\tC\n";
  }
  options.input_path = prefixed_path.string();
  options.source_id_column = "SOURCE_ID";
  options.target_id_column = "TARGET_ID";
  options.target_id_prefix_to_strip = "chr";
  options.output_path = prefixed_index_path.string();
  pgensparsescore::BuildVariantIndex(options);
  {
    pgensparsescore::VariantIndex index(prefixed_index_path.string());
    if (index.Lookup("DRAGEN:chr1:400:G:C") != 0 ||
        index.Lookup("1:400:G:C") != 0 ||
        index.Lookup("chr1:400:G:C").has_value()) {
      throw std::runtime_error("target ID prefix stripping is wrong");
    }
  }

  const auto duplicate_path = directory / "duplicates.tsv";
  const auto duplicate_index = directory / "duplicates.index.bin";
  {
    std::ofstream output(duplicate_path);
    output << "SOURCE_ID\tTARGET_ID\tREF\tALT\n"
           << "duplicate\ttarget-1\tA\tG\n"
           << "duplicate\ttarget-2\tC\tT\n";
  }
  options.input_path = duplicate_path.string();
  options.source_id_column = "SOURCE_ID";
  options.target_id_column = "TARGET_ID";
  options.target_id_prefix_to_strip.clear();
  options.output_path = duplicate_index.string();
  bool rejected = false;
  try {
    pgensparsescore::BuildVariantIndex(options);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  if (!rejected || std::filesystem::exists(duplicate_index) ||
      std::filesystem::exists(duplicate_index.string() + ".tmp")) {
    throw std::runtime_error("duplicate variant-index ID was not rejected");
  }
  std::filesystem::remove_all(directory);
}

void TestScoreFragment() {
  const auto directory = TempPath("-score-fragment");
  std::filesystem::create_directories(directory);
  const auto variants_path = directory / "variants.tsv";
  const auto index_path = directory / "variants.index.bin";
  {
    std::ofstream output(variants_path);
    output << "SOURCE_ID\tTARGET_ID\tREF\tALT\n"
           << "source-v1\ttarget-v1\tA\tG\n"
           << "source-v2\ttarget-v2\tC\tT\n"
           << "source-v3\ttarget-v3\tG\tA\n";
  }
  pgensparsescore::VariantIndexBuildOptions index_options;
  index_options.input_path = variants_path.string();
  index_options.block_size = 2;
  index_options.output_path = index_path.string();
  pgensparsescore::BuildVariantIndex(index_options);

  const auto pvar_path = directory / "reference.pvar";
  const auto frequency_path = directory / "reference.acount";
  const auto support_path = directory / "reference.support.bin";
  {
    std::ofstream output(pvar_path);
    output << "#CHROM\tPOS\tID\tREF\tALT\n"
           << "1\t100\ttarget-v1\tA\tG\n"
           << "1\t200\ttarget-v2\tC\tT\n";
  }
  {
    std::ofstream output(frequency_path);
    output << "#CHROM\tID\tREF\tALT\tREF_CT\tALT_CTS\tOBS_CT\n"
           << "1\ttarget-v1\tA\tG\t5\t1\t6\n";
  }
  pgensparsescore::SupportIndexBuildOptions support_options;
  support_options.variant_index_path = index_path.string();
  support_options.pvar_path = pvar_path.string();
  support_options.frequency_path = frequency_path.string();
  support_options.output_path = support_path.string();
  const auto support_summary =
      pgensparsescore::BuildSupportIndex(support_options);
  if (support_summary.variant_ct != 3 ||
      support_summary.usable_variant_ct != 1 ||
      support_summary.missing_frequency_ct != 1 ||
      support_summary.missing_variant_ct != 1) {
    throw std::runtime_error("support-index summary is wrong");
  }
  {
    pgensparsescore::SupportIndex support(support_path.string());
    if (support.state(0) != pgensparsescore::VariantSupport::kUsable ||
        support.state(1) !=
            pgensparsescore::VariantSupport::kMissingFrequency ||
        support.state(2) !=
            pgensparsescore::VariantSupport::kMissingVariant) {
      throw std::runtime_error("support-index states are wrong");
    }
  }

  const auto weights_a = directory / "weights-a.tsv";
  const auto weights_b = directory / "weights-b.tsv";
  const auto manifest = directory / "manifest.tsv";
  {
    std::ofstream output(weights_a);
    output << "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"
           << "source-v1\tG\tA\t2\n"
           << "source-v2\tC\tT\t3\n"
           << "outside\tA\tG\t5\n"
           << "source-v1\tG\tA\t4\n";
  }
  {
    std::ofstream output(weights_b);
    output << "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"
           << "source-v3\tA\tG\t-1\n"
           << "source-v2\tT\tC\t0\n";
  }
  {
    std::ofstream output(manifest);
    output << "SCORE_INDEX\tSCORE_ID\tCOLUMN_NAME\tPATH\n"
           << "2\tsource:a\tsource__a\tweights-a.tsv\n"
           << "0\tsource:b\tsource__b\tweights-b.tsv\n";
  }
  const auto fragment_path = directory / "scores.fragment.bin";
  pgensparsescore::ScoreFragmentCompileOptions fragment_options;
  fragment_options.manifest_path = manifest.string();
  fragment_options.variant_index_path = index_path.string();
  fragment_options.support_index_path = support_path.string();
  fragment_options.output_path = fragment_path.string();
  const auto summary =
      pgensparsescore::CompileScoreFragment(fragment_options);
  if (summary.score_ct != 2 || summary.catalog_weight_ct != 4 ||
      summary.weight_ct != 2 || summary.missing_variant_weight_ct != 1 ||
      summary.missing_frequency_weight_ct != 1 ||
      summary.input_weight_ct != 6 || summary.zero_weight_ct != 1 ||
      summary.excluded_weight_ct != 1 || summary.duplicate_weight_ct != 1) {
    throw std::runtime_error("score-fragment summary is wrong");
  }
  pgensparsescore::ScoreFragmentReader reader(fragment_path.string());
  if (reader.variant_ct() != 3 || reader.tile_size() != 2 ||
      reader.tile_ct() != 2 || reader.weight_ct() != 2 ||
      reader.scores().size() != 2 ||
      reader.scores()[0].score_id != "source:a" ||
      reader.scores()[1].score_id != "source:b") {
    throw std::runtime_error("score-fragment metadata is wrong");
  }
  const auto tile0 = reader.OpenTile(0);
  if (tile0.variant_ct() != 2 || tile0.referenced_variant_ct() != 1 ||
      tile0.rows().size() != 1 || tile0.rows()[0].local_score_idx() != 0 ||
      tile0.rows()[0].edge_ct() != 2 ||
      tile0.rows()[0].edge(0).local_variant_idx != 0 ||
      tile0.rows()[0].edge(1).local_variant_idx != 0) {
    throw std::runtime_error("score-fragment first tile is wrong");
  }
  ExpectNear(tile0.rows()[0].edge(0).beta_alt +
                 tile0.rows()[0].edge(1).beta_alt,
             6.0,
             "fragment duplicate ALT weights");
  const auto tile1 = reader.OpenTile(1);
  if (tile1.variant_ct() != 1 || tile1.referenced_variant_ct() != 0 ||
      !tile1.rows().empty()) {
    throw std::runtime_error("score-fragment second tile is wrong");
  }
  if (reader.scores()[0].info.catalog_weight_ct != 3 ||
      reader.scores()[0].info.matched_weight_ct != 1 ||
      reader.scores()[0].info.missing_frequency_ct != 1 ||
      reader.scores()[0].info.ref_effect_ct != 1 ||
      reader.scores()[1].info.catalog_weight_ct != 1 ||
      reader.scores()[1].info.missing_variant_ct != 1) {
    throw std::runtime_error("projected fragment QC metadata is wrong");
  }
  std::filesystem::remove_all(directory);
}

void TestMultiFragmentSingleDecodeScoring() {
  constexpr uint32_t sample_ct = 40001;
  constexpr uint32_t variant_ct = 2;
  const auto directory = TempPath("-multi-fragment-score");
  std::filesystem::create_directories(directory);
  const auto pgen_path = directory / "fixture.pgen";
  std::vector<uint64_t> sparse(pgen_rans::PackedWordCt(sample_ct), 0);
  for (uint32_t sample_idx = 0; sample_idx < 4096; ++sample_idx) {
    pgen_rans::SetPackedGenotype(sparse.data(), sample_idx, 1);
  }
  pgen_rans::SetPackedGenotype(sparse.data(), 17, 3);
  std::vector<uint64_t> dense(pgen_rans::PackedWordCt(sample_ct), 0);
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    pgen_rans::SetPackedGenotype(
        dense.data(), sample_idx,
        static_cast<uint8_t>((sample_idx + 1) % 3));
  }
  pgen_rans::SetPackedGenotype(dense.data(), 16, 3);
  pgen_rans::EncodedBlock encoded;
  encoded.first_variant = 0;
  encoded.records.resize(variant_ct);
  std::string error;
  if (!pgen_rans::EncodeSparsePredictorRecord(
          sparse.data(), nullptr, nullptr, sample_ct,
          pgen_rans::RecordMode::kMarginal, 0, 0, &encoded.records[0],
          &error) ||
      !pgen_rans::EncodeRawRecord(dense.data(), sample_ct,
                                 &encoded.records[1], &error)) {
    throw std::runtime_error("cannot encode fragment-scoring fixture: " +
                             error);
  }
  pgen_rans::ContainerWriter writer;
  const pgen_rans::ContainerParams params(
      sample_ct, variant_ct, variant_ct, 1, 4, 12, 1, 1, 2);
  if (!writer.Open(pgen_path.string(), params, {}, &error) ||
      !writer.WriteBlock(encoded, &error) || !writer.Close(&error)) {
    throw std::runtime_error("cannot write fragment-scoring fixture: " +
                             error);
  }

  const auto variants_path = directory / "variants.tsv";
  const auto index_path = directory / "variants.index.bin";
  {
    std::ofstream output(variants_path);
    output << "SOURCE_ID\tTARGET_ID\tREF\tALT\n"
           << "source-v1\ttarget-v1\tA\tG\n"
           << "source-v2\ttarget-v2\tC\tT\n";
  }
  pgensparsescore::VariantIndexBuildOptions index_options;
  index_options.input_path = variants_path.string();
  index_options.block_size = 1;
  index_options.output_path = index_path.string();
  pgensparsescore::BuildVariantIndex(index_options);
  pgensparsescore::VariantIndex index(index_path.string());

  const auto weights_a = directory / "weights-a.tsv";
  const auto weights_b = directory / "weights-b.tsv";
  {
    std::ofstream output(weights_a);
    output << "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"
           << "source-v1\tG\tA\t2\n"
           << "source-v1\tG\tA\t0.125\n"
           << "source-v1\tG\tA\t0.125\n"
           << "source-v1\tG\tA\t0.125\n"
           << "source-v1\tG\tA\t0.125\n"
           << "source-v1\tG\tA\t0.125\n"
           << "source-v1\tG\tA\t0.125\n"
           << "source-v1\tG\tA\t0.125\n"
           << "source-v2\tC\tT\t3\n"
           << "source-v2\tC\tT\t0.5\n";
  }
  {
    std::ofstream output(weights_b);
    output << "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"
           << "source-v1\tA\tG\t4\n"
           << "source-v2\tT\tC\t-1\n";
  }
  const auto combined_manifest = directory / "combined.tsv";
  const auto manifest_a = directory / "fragment-a.tsv";
  const auto manifest_b = directory / "fragment-b.tsv";
  {
    std::ofstream output(combined_manifest);
    output << "SCORE_INDEX\tCOLUMN_NAME\tPATH\n"
           << "0\tscore_a\tweights-a.tsv\n"
           << "1\tscore_b\tweights-b.tsv\n";
  }
  {
    std::ofstream output(manifest_a);
    output << "SCORE_INDEX\tCOLUMN_NAME\tPATH\n"
           << "0\tscore_a\tweights-a.tsv\n";
  }
  {
    std::ofstream output(manifest_b);
    output << "SCORE_INDEX\tCOLUMN_NAME\tPATH\n"
           << "1\tscore_b\tweights-b.tsv\n";
  }
  auto compile = [&](const std::filesystem::path& manifest,
                     const std::filesystem::path& output) {
    pgensparsescore::ScoreFragmentCompileOptions options;
    options.manifest_path = manifest.string();
    options.variant_index_path = index_path.string();
    options.output_path = output.string();
    pgensparsescore::CompileScoreFragment(options);
  };
  const auto combined_path = directory / "combined.fragment.bin";
  const auto fragment_a_path = directory / "a.fragment.bin";
  const auto fragment_b_path = directory / "b.fragment.bin";
  compile(combined_manifest, combined_path);
  compile(manifest_a, fragment_a_path);
  compile(manifest_b, fragment_b_path);

  const std::vector<pgensparsescore::IndexedVariantLocation> locations{
      {0, 0}, {0, 1}};
  pgensparsescore::IndexedFrequencyTable frequencies;
  frequencies.alt_dosage_means = {0.25, 0.75};
  frequencies.matched_row_ct = 2;
  auto score = [&](const std::vector<std::string>& paths,
                   const std::filesystem::path& matrix_path,
                   const std::string& schema_path = "",
                   uint32_t thread_ct = 1) {
    auto loaded =
        pgensparsescore::LoadScoreFragments(paths, index, schema_path);
    pgensparsescore::PgenDosageReader reader(pgen_path.string());
    std::vector<pgensparsescore::PgenDosageReader*> readers{&reader};
    std::vector<double> values;
    pgensparsescore::ScoreRunStats stats;
    {
      pgensparsescore::MappedMatrix matrix(matrix_path.string(), 2, sample_ct);
      stats = pgensparsescore::ScoreFragments(
          index, loaded.fragments, loaded.score_maps, locations, &frequencies,
          pgensparsescore::MissingFrequencyPolicy::kError, readers,
          &loaded.catalog, &matrix, thread_ct);
      values.assign(matrix.Row(0), matrix.Row(0) + 2 * sample_ct);
    }
    std::filesystem::remove(matrix_path);
    return std::make_pair(values, stats);
  };
  const auto combined = score({combined_path.string()},
                              directory / "combined.matrix.bin");
  const auto split = score({fragment_a_path.string(), fragment_b_path.string()},
                           directory / "split.matrix.bin", "", 4);
  if (combined.second.variant_ct != variant_ct ||
      split.second.variant_ct != variant_ct || combined.second.edge_ct != 12 ||
      split.second.edge_ct != 12 ||
      combined.second.sparse_edge_ct + combined.second.dense_edge_ct != 12 ||
      split.second.sparse_edge_ct + split.second.dense_edge_ct != 12 ||
      combined.second.sparse_update_ct + combined.second.dense_update_ct == 0 ||
      split.second.sparse_update_ct + split.second.dense_update_ct == 0 ||
      split.second.score_major_tile_ct != 2 ||
      split.second.score_major_row_ct != 4 ||
      split.second.score_major_maximum_rows_per_tile != 2 ||
      split.second.score_major_maximum_edges_per_tile != 9 ||
      !split.second.score_major_scoring_nanoseconds ||
      split.second.blocked_dense_tile_ct != 1 ||
      split.second.blocked_dense_sample_block_ct != 157 ||
      split.second.blocked_dense_row_ct != 2 ||
      !split.second.blocked_dense_plan_nanoseconds ||
      !split.second.blocked_dense_scoring_nanoseconds ||
      split.second.blocked_dense_sample_block_size != 256 ||
      combined.first != split.first) {
    std::cerr << "fragment stats: dense_tiles="
              << split.second.blocked_dense_tile_ct
              << " dense_blocks="
              << split.second.blocked_dense_sample_block_ct
              << " dense_rows=" << split.second.blocked_dense_row_ct
              << " score_tiles=" << split.second.score_major_tile_ct
              << " score_rows=" << split.second.score_major_row_ct
              << " max_rows="
              << split.second.score_major_maximum_rows_per_tile
              << " max_edges="
              << split.second.score_major_maximum_edges_per_tile
              << " values_equal=" << (combined.first == split.first) << '\n';
    throw std::runtime_error(
        "splitting score fragments changed scores or genotype decode count");
  }
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    double dosage0 = sample_idx < 4096 ? 1.0 : 0.0;
    if (sample_idx == 17) dosage0 = 0.25;
    const double dosage1 = sample_idx == 16
                               ? 0.75
                               : static_cast<double>((sample_idx + 1) % 3);
    ExpectNear(combined.first[sample_idx],
               2.875 * dosage0 + 3.5 * (2.0 - dosage1),
               "multi-fragment score A");
    ExpectNear(combined.first[sample_ct + sample_idx],
               4.0 * (2.0 - dosage0) - dosage1,
               "multi-fragment score B");
  }
  const auto schema_path = directory / "reordered-schema.tsv";
  {
    std::ofstream output(schema_path);
    output << "SCORE_ID\tCOLUMN_NAME\n"
           << "score_b\trenamed_b\n"
           << "score_a\trenamed_a\n";
  }
  const auto reordered = score(
      {fragment_a_path.string(), fragment_b_path.string()},
      directory / "reordered.matrix.bin", schema_path.string());
  for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
    ExpectNear(reordered.first[sample_idx],
               combined.first[sample_ct + sample_idx],
               "stable-ID schema reordered score B");
    ExpectNear(reordered.first[sample_ct + sample_idx],
               combined.first[sample_idx],
               "stable-ID schema reordered score A");
  }
  std::filesystem::remove_all(directory);
}

}  // namespace

int main() {
  try {
    TestDenseKernel();
    TestSparseKernel();
    TestConditionalRansPgenReader();
    TestCatalogOrientation();
    TestRepeatedWeightRowsAreRetained();
    TestVariantMapping();
    TestFilteredPvarKeepsPgenIndexes();
    TestFrequencyParsing();
    TestPfileList();
    TestVariantIndex();
    TestScoreFragment();
    TestMultiFragmentSingleDecodeScoring();
    std::cout << "all unit tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
