// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstdint>
#include <vector>

#include "frequency.h"
#include "mapped_matrix.h"
#include "pgen_reader.h"
#include "types.h"

namespace pgensparsescore {

void ApplyDenseDosage(const double* dosages, uint32_t sample_ct,
                      const std::vector<Edge>& edges, MappedMatrix* matrix);
void ApplySparseDosage(double common, double mean,
                       const uint32_t* sample_ids,
                       const uint16_t* dosage16, uint32_t value_ct,
                       const std::vector<Edge>& edges,
                       std::vector<double>* baselines,
                       MappedMatrix* matrix);
ScoreRunStats ApplyMissingFrequencyPolicy(
    Catalog* catalog, const std::vector<Variant>& variants,
    const FrequencyTable* frequencies, MissingFrequencyPolicy policy);
ScoreRunStats ScoreCatalog(const Catalog& catalog,
                           const std::vector<Variant>& variants,
                           const FrequencyTable* frequencies,
                           MissingFrequencyPolicy missing_frequency_policy,
                           PgenDosageReader* reader, MappedMatrix* matrix);

}  // namespace pgensparsescore
