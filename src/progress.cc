// SPDX-License-Identifier: GPL-3.0-only
#include "progress.h"

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <unistd.h>

#include "io.h"

namespace pgensparsescore {
namespace {

uint64_t PeakRssBytes() {
  struct rusage usage {};
  if (getrusage(RUSAGE_SELF, &usage)) {
    return 0;
  }
#ifdef __APPLE__
  return static_cast<uint64_t>(usage.ru_maxrss);
#else
  return static_cast<uint64_t>(usage.ru_maxrss) * 1024;
#endif
}

uint64_t CurrentRssBytes() {
#ifdef __linux__
  std::ifstream input("/proc/self/statm");
  uint64_t total_pages = 0;
  uint64_t resident_pages = 0;
  if (input >> total_pages >> resident_pages) {
    static_cast<void>(total_pages);
    const long page_size = sysconf(_SC_PAGESIZE);
    if (page_size > 0) {
      return resident_pages * static_cast<uint64_t>(page_size);
    }
  }
#endif
  return 0;
}

}  // namespace

ProgressReporter::ProgressReporter(const std::string& path,
                                   uint32_t interval_seconds)
    : output_(path, std::ios::out | std::ios::trunc),
      started_(std::chrono::steady_clock::now()),
      interval_(interval_seconds),
      sequence_(0) {
  if (!output_) {
    throw std::runtime_error("cannot create progress log " + path);
  }
  next_event_ = started_;
}

void ProgressReporter::Event(const std::string& operation,
                             const std::string& phase,
                             const ProgressCounters& counters,
                             const ProgressDetails& details) {
  if (!enabled()) {
    return;
  }
  Write(operation, phase, counters, details);
}

void ProgressReporter::MaybeEvent(const std::string& operation,
                                  const std::string& phase,
                                  const ProgressCounters& counters,
                                  const ProgressDetails& details) {
  if (!enabled() || std::chrono::steady_clock::now() < next_event_) {
    return;
  }
  Write(operation, phase, counters, details);
}

void ProgressReporter::Write(const std::string& operation,
                             const std::string& phase,
                             const ProgressCounters& counters,
                             const ProgressDetails& details) {
  const auto now_steady = std::chrono::steady_clock::now();
  const auto now_system = std::chrono::system_clock::now();
  const uint64_t timestamp_unix_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
          now_system.time_since_epoch())
          .count());
  const uint64_t elapsed_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(now_steady -
                                                            started_)
          .count());
  next_event_ = now_steady + interval_;

  output_ << "{\"schema_version\":1,\"sequence\":" << sequence_++
          << ",\"timestamp_unix_ms\":" << timestamp_unix_ms
          << ",\"elapsed_ms\":" << elapsed_ms << ",\"operation\":\""
          << JsonEscape(operation) << "\",\"phase\":\""
          << JsonEscape(phase) << "\",\"pid\":" << getpid()
          << ",\"rss_bytes\":" << CurrentRssBytes()
          << ",\"peak_rss_bytes\":" << PeakRssBytes();
  for (const auto& [name, value] : counters) {
    output_ << ",\"" << JsonEscape(name) << "\":" << value;
  }
  for (const auto& [name, value] : details) {
    output_ << ",\"" << JsonEscape(name) << "\":\""
            << JsonEscape(value) << "\"";
  }
  output_ << "}\n";
  output_.flush();
  if (!output_) {
    throw std::runtime_error("cannot write progress log");
  }

  std::cerr << "progress " << operation << ':' << phase;
  for (const auto& [name, value] : counters) {
    std::cerr << ' ' << name << '=' << value;
  }
  std::cerr << '\n';
}

}  // namespace pgensparsescore
