// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <chrono>
#include <cstdint>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace pgensparsescore {

using ProgressCounters = std::vector<std::pair<std::string, uint64_t>>;
using ProgressDetails = std::vector<std::pair<std::string, std::string>>;

class ProgressReporter {
 public:
  ProgressReporter() = default;
  ProgressReporter(const std::string& path, uint32_t interval_seconds);

  bool enabled() const { return output_.is_open(); }
  void Event(const std::string& operation, const std::string& phase,
             const ProgressCounters& counters = {},
             const ProgressDetails& details = {});
  void MaybeEvent(const std::string& operation, const std::string& phase,
                  const ProgressCounters& counters = {},
                  const ProgressDetails& details = {});

 private:
  void Write(const std::string& operation, const std::string& phase,
             const ProgressCounters& counters,
             const ProgressDetails& details);

  std::ofstream output_;
  std::chrono::steady_clock::time_point started_;
  std::chrono::steady_clock::time_point next_event_;
  std::chrono::seconds interval_{30};
  uint64_t sequence_ = 0;
};

}  // namespace pgensparsescore
