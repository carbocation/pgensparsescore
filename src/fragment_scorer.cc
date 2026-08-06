// SPDX-License-Identifier: GPL-3.0-only
#include "fragment_scorer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <optional>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "io.h"
#include "scorer.h"

namespace pgensparsescore {
namespace {

using Header = std::unordered_map<std::string, size_t>;

Header MakeHeader(const std::vector<std::string>& fields,
                  const std::string& path) {
  Header result;
  for (size_t idx = 0; idx < fields.size(); ++idx) {
    std::string name = fields[idx];
    if (!name.empty() && name.front() == '#') name.erase(0, 1);
    if (!result.emplace(name, idx).second) {
      throw std::runtime_error(path + " has duplicate column " + name);
    }
  }
  return result;
}

size_t RequireColumn(const Header& header, const std::string& name,
                     const std::string& path) {
  const auto iter = header.find(name);
  if (iter == header.end()) {
    throw std::runtime_error(path + " is missing column " + name);
  }
  return iter->second;
}

struct QueueEntry {
  uint32_t ordinal;
  uint32_t fragment_idx;
};

struct LaterOrdinal {
  bool operator()(const QueueEntry& lhs, const QueueEntry& rhs) const {
    if (lhs.ordinal != rhs.ordinal) return lhs.ordinal > rhs.ordinal;
    return lhs.fragment_idx > rhs.fragment_idx;
  }
};

}  // namespace

std::vector<std::string> ReadScoreFragmentList(const std::string& path) {
  LineReader reader(path);
  std::string line;
  if (!reader.GetLine(&line)) throw std::runtime_error(path + " is empty");
  const Header header = MakeHeader(SplitTabs(line), path);
  size_t fragment_idx = static_cast<size_t>(-1);
  for (const char* name : {"FRAGMENT", "PATH"}) {
    const auto iter = header.find(name);
    if (iter != header.end()) {
      fragment_idx = iter->second;
      break;
    }
  }
  if (fragment_idx == static_cast<size_t>(-1)) {
    throw std::runtime_error(path + " is missing FRAGMENT or PATH column");
  }
  const auto base = std::filesystem::absolute(path).parent_path();
  std::vector<std::string> result;
  std::unordered_set<std::string> seen;
  uint64_t line_number = 1;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty()) continue;
    const auto fields = SplitTabs(line);
    if (fields.size() <= fragment_idx || fields[fragment_idx].empty()) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has no fragment path");
    }
    std::filesystem::path resolved(fields[fragment_idx]);
    if (resolved.is_relative()) resolved = base / resolved;
    const std::string normalized = resolved.lexically_normal().string();
    if (!seen.insert(normalized).second) {
      throw std::runtime_error(path + " contains duplicate fragment " +
                               fields[fragment_idx]);
    }
    result.push_back(normalized);
  }
  if (result.empty()) throw std::runtime_error(path + " has no fragments");
  return result;
}

IndexedPvarStats AddIndexedPvar(
    const std::string& path, uint32_t input_idx, const VariantIndex& index,
    std::vector<IndexedVariantLocation>* locations,
    ProgressReporter* progress) {
  if (locations->size() != index.variant_ct()) {
    throw std::runtime_error("variant-location vector does not match index");
  }
  LineReader reader(path);
  std::string line;
  Header header;
  size_t id_idx = 0;
  size_t ref_idx = 0;
  size_t alt_idx = 0;
  size_t maximum = 0;
  uint64_t line_number = 0;
  IndexedPvarStats stats;
  while (reader.GetLine(&line)) {
    ++line_number;
    if (line.empty() || line.rfind("##", 0) == 0) continue;
    if (header.empty()) {
      header = MakeHeader(SplitTabs(line), path);
      id_idx = RequireColumn(header, "ID", path);
      ref_idx = RequireColumn(header, "REF", path);
      alt_idx = RequireColumn(header, "ALT", path);
      maximum = std::max({id_idx, ref_idx, alt_idx});
      continue;
    }
    if (stats.row_ct == UINT32_MAX) {
      throw std::runtime_error(path + " exceeds the supported PVAR row count");
    }
    const uint32_t pgen_variant_idx = stats.row_ct++;
    const auto fields = SplitTabs(line);
    if (fields.size() <= maximum) {
      throw std::runtime_error(path + ": line " +
                               std::to_string(line_number) +
                               " has too few fields");
    }
    const auto ordinal = index.Lookup(fields[id_idx]);
    if (!ordinal) {
      if (progress && !(stats.row_ct % 1000000)) {
        progress->MaybeEvent("score", "read_pvar",
                             {{"pgen_input_index", input_idx},
                              {"pvar_rows_scanned", stats.row_ct},
                              {"indexed_variants_matched", stats.matched_variant_ct}},
                             {{"pvar", path}});
      }
      continue;
    }
    if (fields[alt_idx].find(',') != std::string::npos ||
        fields[ref_idx] != index.ref(*ordinal) ||
        fields[alt_idx] != index.alt(*ordinal)) {
      throw std::runtime_error(
          path + ": line " + std::to_string(line_number) +
          " alleles disagree with the variant index for " + fields[id_idx]);
    }
    auto& location = locations->at(*ordinal);
    if (location.present()) {
      throw std::runtime_error(
          "indexed variant occurs in more than one PVAR row: " +
          fields[id_idx]);
    }
    location.input_idx = input_idx;
    location.pgen_variant_idx = pgen_variant_idx;
    ++stats.matched_variant_ct;
    if (progress && !(stats.row_ct % 1000000)) {
      progress->MaybeEvent("score", "read_pvar",
                           {{"pgen_input_index", input_idx},
                            {"pvar_rows_scanned", stats.row_ct},
                            {"indexed_variants_matched", stats.matched_variant_ct}},
                           {{"pvar", path}});
    }
  }
  if (header.empty()) throw std::runtime_error(path + " has no header");
  if (!stats.row_ct) throw std::runtime_error(path + " has no variants");
  if (progress) {
    progress->Event("score", "pvar_loaded",
                    {{"pgen_input_index", input_idx},
                     {"pvar_rows", stats.row_ct},
                     {"indexed_variants_matched", stats.matched_variant_ct}},
                    {{"pvar", path}});
  }
  return stats;
}

LoadedScoreFragments LoadScoreFragments(
    const std::vector<std::string>& paths, const VariantIndex& index,
    const std::string& score_schema_path, ProgressReporter* progress) {
  if (paths.empty()) throw std::runtime_error("no score fragments supplied");
  LoadedScoreFragments result;
  result.fragments.reserve(paths.size());
  uint64_t score_ct = 0;
  for (const auto& path : paths) {
    result.fragments.emplace_back(path);
    const auto& fragment = result.fragments.back();
    if (fragment.variant_ct() != index.variant_ct() ||
        fragment.block_size() != index.block_size() ||
        fragment.block_ct() != index.block_ct() ||
        fragment.signature_lo() != index.signature_lo() ||
        fragment.signature_hi() != index.signature_hi()) {
      throw std::runtime_error(path +
                               " was built for a different variant index");
    }
    score_ct += fragment.scores().size();
    result.weight_ct += fragment.weight_ct();
    result.file_byte_ct += fragment.file_bytes();
    if (progress) {
      progress->Event("score", "fragment_loaded",
                      {{"fragments_loaded", result.fragments.size()},
                       {"fragments_total", paths.size()},
                       {"scores_loaded", score_ct},
                       {"weights_available", result.weight_ct},
                       {"fragment_bytes", result.file_byte_ct}},
                      {{"fragment", path}});
    }
  }
  if (!score_ct || score_ct > UINT32_MAX) {
    throw std::runtime_error("score fragments contain too many scores");
  }

  std::vector<std::string> schema_ids;
  std::vector<std::string> schema_columns;
  if (!score_schema_path.empty()) {
    LineReader reader(score_schema_path);
    std::string line;
    if (!reader.GetLine(&line)) {
      throw std::runtime_error(score_schema_path + " is empty");
    }
    const Header header = MakeHeader(SplitTabs(line), score_schema_path);
    size_t id_idx = static_cast<size_t>(-1);
    for (const char* name : {"SCORE_ID", "SCORE", "COLUMN_NAME"}) {
      const auto iter = header.find(name);
      if (iter != header.end()) {
        id_idx = iter->second;
        break;
      }
    }
    if (id_idx == static_cast<size_t>(-1)) {
      throw std::runtime_error(
          score_schema_path + " is missing SCORE_ID, SCORE, or COLUMN_NAME");
    }
    const auto column_iter = header.find("COLUMN_NAME");
    const size_t column_idx =
        column_iter == header.end() ? id_idx : column_iter->second;
    const size_t maximum = std::max(id_idx, column_idx);
    std::unordered_set<std::string> seen_schema_ids;
    std::unordered_set<std::string> seen_schema_columns;
    uint64_t line_number = 1;
    while (reader.GetLine(&line)) {
      ++line_number;
      if (line.empty()) continue;
      const auto fields = SplitTabs(line);
      if (fields.size() <= maximum || fields[id_idx].empty() ||
          fields[column_idx].empty() ||
          !seen_schema_ids.insert(fields[id_idx]).second ||
          !seen_schema_columns.insert(fields[column_idx]).second) {
        throw std::runtime_error(
            score_schema_path + ": line " + std::to_string(line_number) +
            " has an empty or duplicate score ID/column");
      }
      schema_ids.push_back(fields[id_idx]);
      schema_columns.push_back(fields[column_idx]);
    }
    if (schema_ids.size() != score_ct) {
      throw std::runtime_error(
          "score schema and fragments contain different score counts");
    }
  }

  std::unordered_map<std::string, uint32_t> output_by_id;
  if (!schema_ids.empty()) {
    output_by_id.reserve(schema_ids.size());
    for (uint32_t idx = 0; idx < schema_ids.size(); ++idx) {
      output_by_id.emplace(schema_ids[idx], idx);
    }
  }
  result.catalog.scores.resize(static_cast<size_t>(score_ct));
  result.catalog.intercepts.assign(static_cast<size_t>(score_ct), 0.0);
  result.score_maps.resize(result.fragments.size());
  std::vector<bool> output_seen(score_ct, false);
  std::unordered_set<std::string> fragment_ids;
  std::unordered_set<std::string> output_columns;
  uint32_t next_output_idx = 0;
  for (uint32_t fragment_idx = 0; fragment_idx < result.fragments.size();
       ++fragment_idx) {
    const auto& fragment = result.fragments[fragment_idx];
    auto& score_map = result.score_maps[fragment_idx];
    score_map.reserve(fragment.scores().size());
    for (const auto& fragment_score : fragment.scores()) {
      if (!fragment_ids.insert(fragment_score.score_id).second) {
        throw std::runtime_error(
            "score fragments contain duplicate stable score ID " +
            fragment_score.score_id);
      }
      uint32_t output_idx = next_output_idx++;
      std::string output_column = fragment_score.info.id;
      if (!schema_ids.empty()) {
        const auto schema = output_by_id.find(fragment_score.score_id);
        if (schema == output_by_id.end()) {
          throw std::runtime_error("score schema has no row for fragment score " +
                                   fragment_score.score_id);
        }
        output_idx = schema->second;
        output_column = schema_columns[output_idx];
      }
      if (output_seen[output_idx] ||
          !output_columns.insert(output_column).second) {
        throw std::runtime_error(
            "score fragments resolve to a duplicate output column");
      }
      output_seen[output_idx] = true;
      score_map.push_back(output_idx);
      result.catalog.scores[output_idx] = fragment_score.info;
      auto& score = result.catalog.scores[output_idx];
      score.id = output_column;
      score.matched_weight_ct = 0;
      score.missing_variant_ct = 0;
      score.missing_frequency_ct = 0;
      score.ref_effect_ct = 0;
      score.alt_effect_ct = 0;
      score.ref_effect_intercept = 0.0;
    }
  }
  if (std::find(output_seen.begin(), output_seen.end(), false) !=
      output_seen.end()) {
    throw std::runtime_error("score schema contains a score absent from fragments");
  }
  if (progress) {
    progress->Event("score", "fragments_ready",
                    {{"fragments", result.fragments.size()},
                     {"scores", result.catalog.scores.size()},
                     {"weights_available", result.weight_ct},
                     {"fragment_bytes", result.file_byte_ct}});
  }
  return result;
}

ScoreRunStats ScoreFragments(
    const VariantIndex& index,
    const std::vector<ScoreFragmentReader>& fragments,
    const std::vector<std::vector<uint32_t>>& score_maps,
    const std::vector<IndexedVariantLocation>& locations,
    const IndexedFrequencyTable* frequencies,
    MissingFrequencyPolicy missing_frequency_policy,
    const std::vector<PgenDosageReader*>& readers, Catalog* catalog,
    MappedMatrix* matrix, ProgressReporter* progress) {
  if (locations.size() != index.variant_ct() || fragments.empty() ||
      score_maps.size() != fragments.size() ||
      readers.empty() || matrix->row_ct() != catalog->scores.size() ||
      matrix->column_ct() != readers.front()->sample_ct()) {
    throw std::runtime_error("fragment scoring inputs have incompatible shapes");
  }
  for (const auto* reader : readers) {
    if (!reader || reader->sample_ct() != matrix->column_ct()) {
      throw std::runtime_error("fragment PGEN readers have different samples");
    }
  }
  if (frequencies && frequencies->alt_dosage_means.size() != index.variant_ct()) {
    throw std::runtime_error("indexed frequencies do not match variant index");
  }
  if (!frequencies && missing_frequency_policy != MissingFrequencyPolicy::kCohort) {
    throw std::runtime_error(
        "error and omit missing-frequency policies require frequencies");
  }

  ScoreRunStats stats;
  std::vector<double> baselines(catalog->scores.size(), 0.0);
  std::vector<Edge> combined_edges;
  uint64_t variant_groups_processed = 0;

  for (uint32_t block_idx = 0; block_idx < index.block_ct(); ++block_idx) {
    std::vector<ScoreFragmentBlockCursor> cursors;
    cursors.reserve(fragments.size());
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, LaterOrdinal> queue;
    for (uint32_t fragment_idx = 0; fragment_idx < fragments.size();
         ++fragment_idx) {
      cursors.push_back(fragments[fragment_idx].OpenBlock(block_idx));
      if (!cursors.back().done()) {
        queue.push({cursors.back().ordinal(), fragment_idx});
      }
    }
    while (!queue.empty()) {
      const uint32_t ordinal = queue.top().ordinal;
      combined_edges.clear();
      while (!queue.empty() && queue.top().ordinal == ordinal) {
        const QueueEntry entry = queue.top();
        queue.pop();
        auto& cursor = cursors[entry.fragment_idx];
        const size_t edge_begin = combined_edges.size();
        cursor.AppendEdges(&combined_edges);
        const auto& score_map = score_maps[entry.fragment_idx];
        for (size_t edge_idx = edge_begin; edge_idx < combined_edges.size();
             ++edge_idx) {
          const uint32_t local_score_idx = combined_edges[edge_idx].score_idx;
          if (local_score_idx >= score_map.size()) {
            throw std::runtime_error(
                "fragment edge has invalid local score index");
          }
          combined_edges[edge_idx].score_idx = score_map[local_score_idx];
        }
        cursor.Next();
        if (!cursor.done()) {
          queue.push({cursor.ordinal(), entry.fragment_idx});
        }
      }
      ++variant_groups_processed;
      const auto& location = locations[ordinal];
      if (!location.present()) {
        for (const auto& edge : combined_edges) {
          if (edge.score_idx >= catalog->scores.size()) {
            throw std::runtime_error("fragment edge has invalid SCORE_INDEX");
          }
          ++catalog->scores[edge.score_idx].missing_variant_ct;
        }
        continue;
      }
      for (const auto& edge : combined_edges) {
        if (edge.score_idx >= catalog->scores.size()) {
          throw std::runtime_error("fragment edge has invalid SCORE_INDEX");
        }
        auto& score = catalog->scores[edge.score_idx];
        ++score.matched_weight_ct;
        if (edge.ref_effect) {
          ++score.ref_effect_ct;
        } else {
          ++score.alt_effect_ct;
        }
      }

      std::optional<double> imputation_mean;
      bool omit = false;
      if (frequencies) {
        const double mean = frequencies->alt_dosage_means[ordinal];
        if (std::isnan(mean)) {
          ++stats.missing_frequency_variant_ct;
          if (missing_frequency_policy == MissingFrequencyPolicy::kError) {
            throw std::runtime_error(
                "frequency file has no usable row for indexed variant ordinal " +
                std::to_string(ordinal));
          }
          if (missing_frequency_policy == MissingFrequencyPolicy::kOmit) {
            omit = true;
            ++stats.omitted_frequency_variant_ct;
            stats.omitted_frequency_edge_ct += combined_edges.size();
            for (const auto& edge : combined_edges) {
              ++catalog->scores[edge.score_idx].missing_frequency_ct;
            }
          } else {
            ++stats.cohort_frequency_variant_ct;
          }
        } else {
          imputation_mean = mean;
          ++stats.external_frequency_variant_ct;
        }
      } else {
        ++stats.cohort_frequency_variant_ct;
      }
      if (omit) continue;

      for (const auto& edge : combined_edges) {
        if (edge.ref_effect) {
          const double intercept = -2.0 * edge.beta_alt;
          baselines[edge.score_idx] += intercept;
          catalog->intercepts[edge.score_idx] += intercept;
          catalog->scores[edge.score_idx].ref_effect_intercept += intercept;
        }
      }
      if (location.input_idx >= readers.size()) {
        throw std::runtime_error("indexed PVAR location has invalid input index");
      }
      PgenDosageReader* reader = readers[location.input_idx];
      const DosageView dosage =
          reader->Read(location.pgen_variant_idx, imputation_mean);
      ++stats.variant_ct;
      stats.edge_ct += combined_edges.size();
      stats.imputed_value_ct += dosage.missing_ct;
      if (dosage.sparse) {
        ++stats.sparse_variant_ct;
        ApplySparseDosage(dosage.common, dosage.mean, dosage.sparse_sample_ids,
                          dosage.sparse_dosage16, dosage.sparse_value_ct,
                          combined_edges, &baselines, matrix);
      } else {
        ++stats.dense_variant_ct;
        ApplyDenseDosage(dosage.dense_values, reader->sample_ct(),
                         combined_edges, matrix);
      }
      if (progress && !(stats.variant_ct % 10000)) {
        progress->MaybeEvent(
            "score", "score_fragment_blocks",
            {{"blocks_processed", block_idx},
             {"blocks_total", index.block_ct()},
             {"variant_groups_processed", variant_groups_processed},
             {"genotype_decodes", stats.variant_ct},
             {"weight_edges_processed", stats.edge_ct},
             {"sparse_variants", stats.sparse_variant_ct},
             {"dense_variants", stats.dense_variant_ct},
             {"imputed_values", stats.imputed_value_ct}});
      }
      if (!(stats.variant_ct % 100000)) {
        std::cerr << "decoded " << stats.variant_ct
                  << " indexed variants across " << fragments.size()
                  << " score fragments\n";
      }
    }
    if (progress) {
      progress->MaybeEvent(
          "score", "score_fragment_blocks",
          {{"blocks_processed", block_idx + 1},
           {"blocks_total", index.block_ct()},
           {"variant_groups_processed", variant_groups_processed},
           {"genotype_decodes", stats.variant_ct},
           {"weight_edges_processed", stats.edge_ct},
           {"sparse_variants", stats.sparse_variant_ct},
           {"dense_variants", stats.dense_variant_ct},
           {"imputed_values", stats.imputed_value_ct}});
    }
  }
  for (uint32_t score_idx = 0; score_idx < catalog->scores.size(); ++score_idx) {
    double* row = matrix->Row(score_idx);
    const double baseline = baselines[score_idx];
    for (uint32_t sample_idx = 0; sample_idx < matrix->column_ct(); ++sample_idx) {
      row[sample_idx] += baseline;
    }
  }
  matrix->Flush();
  if (progress) {
    progress->Event(
        "score", "fragment_scoring_complete",
        {{"fragments", fragments.size()},
         {"blocks_processed", index.block_ct()},
         {"variant_groups_processed", variant_groups_processed},
         {"genotype_decodes", stats.variant_ct},
         {"weight_edges_processed", stats.edge_ct},
         {"sparse_variants", stats.sparse_variant_ct},
         {"dense_variants", stats.dense_variant_ct},
         {"imputed_values", stats.imputed_value_ct}});
  }
  return stats;
}

}  // namespace pgensparsescore
