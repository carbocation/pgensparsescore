// SPDX-License-Identifier: GPL-3.0-only
#include "scorer.h"

#include <iostream>
#include <optional>
#include <stdexcept>

namespace pgensparsescore {

void ApplyDenseDosage(const double* dosages, uint32_t sample_ct,
                      const std::vector<Edge>& edges, MappedMatrix* matrix) {
  for (const Edge& edge : edges) {
    double* row = matrix->Row(edge.score_idx);
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      row[sample_idx] += edge.beta_alt * dosages[sample_idx];
    }
  }
}

void ApplySparseDosage(double common, double mean,
                       const uint32_t* sample_ids,
                       const uint16_t* dosage16, uint32_t value_ct,
                       const std::vector<Edge>& edges,
                       std::vector<double>* baselines,
                       MappedMatrix* matrix) {
  for (const Edge& edge : edges) {
    (*baselines)[edge.score_idx] += edge.beta_alt * common;
    double* row = matrix->Row(edge.score_idx);
    for (uint32_t value_idx = 0; value_idx < value_ct; ++value_idx) {
      const double dosage = dosage16[value_idx] == UINT16_MAX
                                ? mean
                                : static_cast<double>(dosage16[value_idx]) /
                                      16384.0;
      row[sample_ids[value_idx]] += edge.beta_alt * (dosage - common);
    }
  }
}

ScoreRunStats ScoreCatalog(const Catalog& catalog,
                           const std::vector<Variant>& variants,
                           const FrequencyTable* frequencies,
                           bool error_on_missing_frequency,
                           PgenDosageReader* reader, MappedMatrix* matrix) {
  if (matrix->row_ct() != catalog.scores.size() ||
      matrix->column_ct() != reader->sample_ct()) {
    throw std::runtime_error("score matrix shape does not match catalog/PGEN");
  }
  std::vector<double> baselines = catalog.intercepts;
  ScoreRunStats stats;
  stats.variant_ct = catalog.variants.size();
  uint64_t processed = 0;
  for (const VariantEdges& variant : catalog.variants) {
    const Variant& variant_metadata = variants.at(variant.variant_idx);
    std::optional<double> imputation_mean;
    if (frequencies) {
      const auto frequency = frequencies->find(variant_metadata.id);
      if (frequency == frequencies->end()) {
        ++stats.missing_frequency_variant_ct;
        if (error_on_missing_frequency) {
          throw std::runtime_error("frequency file has no row for scored variant " +
                                   variant_metadata.id);
        }
        ++stats.cohort_frequency_variant_ct;
      } else {
        if (frequency->second.ref != variant_metadata.ref ||
            frequency->second.alt != variant_metadata.alt) {
          throw std::runtime_error(
              "frequency alleles disagree with PVAR for " +
              variant_metadata.id + " (PVAR " + variant_metadata.ref + "/" +
              variant_metadata.alt + ", frequency " + frequency->second.ref +
              "/" + frequency->second.alt + ")");
        }
        imputation_mean = frequency->second.alt_dosage_mean;
        ++stats.external_frequency_variant_ct;
      }
    } else {
      ++stats.cohort_frequency_variant_ct;
    }
    const DosageView dosage =
        reader->Read(variant_metadata.pgen_variant_idx, imputation_mean);
    stats.edge_ct += variant.edges.size();
    stats.imputed_value_ct += dosage.missing_ct;
    if (dosage.sparse) {
      ++stats.sparse_variant_ct;
      ApplySparseDosage(dosage.common, dosage.mean,
                        dosage.sparse_sample_ids, dosage.sparse_dosage16,
                        dosage.sparse_value_ct, variant.edges, &baselines,
                        matrix);
    } else {
      ++stats.dense_variant_ct;
      ApplyDenseDosage(dosage.dense_values, reader->sample_ct(), variant.edges,
                       matrix);
    }
    ++processed;
    if (processed % 100000 == 0) {
      std::cerr << "processed " << processed << "/" << stats.variant_ct
                << " scored variants\n";
    }
  }
  for (uint32_t score_idx = 0; score_idx < catalog.scores.size(); ++score_idx) {
    double* row = matrix->Row(score_idx);
    const double baseline = baselines[score_idx];
    for (uint32_t sample_idx = 0; sample_idx < reader->sample_ct();
         ++sample_idx) {
      row[sample_idx] += baseline;
    }
  }
  matrix->Flush();
  return stats;
}

}  // namespace pgensparsescore
