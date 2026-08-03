// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include <zlib.h>
#include <zstd.h>

namespace pgensparsescore {

class LineReader {
 public:
  explicit LineReader(const std::string& path);
  ~LineReader();

  LineReader(const LineReader&) = delete;
  LineReader& operator=(const LineReader&) = delete;

  bool GetLine(std::string* line);

 private:
  enum class Mode { kPlain, kGzip, kZstd };

  bool GetZstdLine(std::string* line);

  Mode mode_ = Mode::kPlain;
  std::ifstream plain_;
  gzFile gz_ = nullptr;
  ZSTD_DStream* zstd_ = nullptr;
  std::vector<char> zstd_input_;
  std::vector<char> zstd_output_;
  size_t zstd_input_pos_ = 0;
  size_t zstd_input_size_ = 0;
  size_t zstd_hint_ = 1;
  bool zstd_input_eof_ = false;
  bool zstd_finished_ = false;
  std::string zstd_pending_;
  size_t zstd_pending_pos_ = 0;
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
