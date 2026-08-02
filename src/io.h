// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include <zlib.h>

namespace pgensparsescore {

class LineReader {
 public:
  explicit LineReader(const std::string& path);
  ~LineReader();

  LineReader(const LineReader&) = delete;
  LineReader& operator=(const LineReader&) = delete;

  bool GetLine(std::string* line);

 private:
  bool gzip_ = false;
  std::ifstream plain_;
  gzFile gz_ = nullptr;
  std::string path_;
};

class GzipWriter {
 public:
  explicit GzipWriter(const std::string& path);
  ~GzipWriter();

  GzipWriter(const GzipWriter&) = delete;
  GzipWriter& operator=(const GzipWriter&) = delete;

  void Write(std::string_view value);
  void Close();

 private:
  gzFile gz_ = nullptr;
  std::string path_;
};

std::vector<std::string> SplitTabs(const std::string& line);
std::string JsonEscape(const std::string& value);

}  // namespace pgensparsescore
