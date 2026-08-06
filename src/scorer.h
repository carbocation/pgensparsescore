// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <vector>

#include "frequency.h"
#include "mapped_matrix.h"
#include "pgen_reader.h"
#include "progress.h"
#include "types.h"

namespace pgensparsescore {

void ApplyDenseDosage(const double* dosages, uint32_t sample_ct,
                      const std::vector<Edge>& edges, MappedMatrix* matrix);
void ApplyDenseDosageRange(const double* dosages, uint32_t sample_ct,
                           const std::vector<Edge>& edges, size_t edge_begin,
                           size_t edge_end, MappedMatrix* matrix);
void ApplyDenseDosagePartition(const double* dosages, uint32_t sample_ct,
                               const std::vector<Edge>& edges,
                               uint32_t partition_idx,
                               uint32_t partition_ct, MappedMatrix* matrix);
void ApplySparseDosage(double common, double mean,
                       const uint32_t* sample_ids,
                       const uint16_t* dosage16, uint32_t value_ct,
                       const std::vector<Edge>& edges,
                       std::vector<double>* baselines,
                       MappedMatrix* matrix);
void ApplySparseDosageRange(double common, double mean,
                            const uint32_t* sample_ids,
                            const uint16_t* dosage16, uint32_t value_ct,
                            const std::vector<Edge>& edges, size_t edge_begin,
                            size_t edge_end, std::vector<double>* baselines,
                            MappedMatrix* matrix);
void ApplySparseDosagePartition(double common, double mean,
                                const uint32_t* sample_ids,
                                const uint16_t* dosage16, uint32_t value_ct,
                                const std::vector<Edge>& edges,
                                uint32_t partition_idx,
                                uint32_t partition_ct,
                                std::vector<double>* baselines,
                                MappedMatrix* matrix);
ScoreRunStats ApplyMissingFrequencyPolicy(
    Catalog* catalog, const std::vector<Variant>& variants,
    const FrequencyTable* frequencies, MissingFrequencyPolicy policy);
ScoreRunStats ScoreCatalog(const Catalog& catalog,
                           const std::vector<Variant>& variants,
                           const FrequencyTable* frequencies,
                           MissingFrequencyPolicy missing_frequency_policy,
                           PgenDosageReader* reader, MappedMatrix* matrix,
                           ProgressReporter* progress = nullptr);

}  // namespace pgensparsescore
