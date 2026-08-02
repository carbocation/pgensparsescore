// SPDX-License-Identifier: GPL-3.0-only
#include "io.h"

#include <array>
#include <limits>
#include <stdexcept>

namespace pgensparsescore {

namespace {

bool EndsWith(const std::string& value, const std::string& suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) ==
             0;
}

}  // namespace

LineReader::LineReader(const std::string& path) : path_(path) {
  gzip_ = EndsWith(path, ".gz");
  if (gzip_) {
    gz_ = gzopen(path.c_str(), "rb");
    if (!gz_) {
      throw std::runtime_error("cannot open " + path);
    }
  } else {
    plain_.open(path);
    if (!plain_) {
      throw std::runtime_error("cannot open " + path);
    }
  }
}

LineReader::~LineReader() {
  if (gz_) {
    gzclose(gz_);
  }
}

bool LineReader::GetLine(std::string* line) {
  line->clear();
  if (!gzip_) {
    if (!std::getline(plain_, *line)) {
      return false;
    }
  } else {
    std::array<char, 65536> buffer{};
    while (true) {
      char* result = gzgets(gz_, buffer.data(), static_cast<int>(buffer.size()));
      if (!result) {
        if (line->empty()) {
          int error_number = Z_OK;
          const char* message = gzerror(gz_, &error_number);
          if (error_number != Z_OK && error_number != Z_STREAM_END) {
            throw std::runtime_error("error reading " + path_ + ": " +
                                     (message ? message : "unknown zlib error"));
          }
          return false;
        }
        break;
      }
      line->append(buffer.data());
      if (!line->empty() && line->back() == '\n') {
        break;
      }
    }
    if (!line->empty() && line->back() == '\n') {
      line->pop_back();
    }
  }
  if (!line->empty() && line->back() == '\r') {
    line->pop_back();
  }
  return true;
}

GzipWriter::GzipWriter(const std::string& path) : path_(path) {
  gz_ = gzopen(path.c_str(), "wb6");
  if (!gz_) {
    throw std::runtime_error("cannot create " + path);
  }
}

GzipWriter::~GzipWriter() {
  if (gz_) {
    gzclose(gz_);
  }
}

void GzipWriter::Write(std::string_view value) {
  size_t offset = 0;
  while (offset < value.size()) {
    const size_t chunk_size = std::min(
        value.size() - offset,
        static_cast<size_t>(std::numeric_limits<int>::max()));
    const int written = gzwrite(gz_, value.data() + offset,
                                static_cast<unsigned int>(chunk_size));
    if (written <= 0 || static_cast<size_t>(written) != chunk_size) {
      int error_number = Z_OK;
      const char* message = gzerror(gz_, &error_number);
      throw std::runtime_error("error writing " + path_ + ": " +
                               (message ? message : "unknown zlib error"));
    }
    offset += chunk_size;
  }
}

void GzipWriter::Close() {
  if (!gz_) {
    return;
  }
  const int result = gzclose(gz_);
  gz_ = nullptr;
  if (result != Z_OK) {
    throw std::runtime_error("error closing " + path_);
  }
}

std::vector<std::string> SplitTabs(const std::string& line) {
  std::vector<std::string> fields;
  size_t begin = 0;
  while (true) {
    const size_t end = line.find('\t', begin);
    if (end == std::string::npos) {
      fields.emplace_back(line.substr(begin));
      return fields;
    }
    fields.emplace_back(line.substr(begin, end - begin));
    begin = end + 1;
  }
}

std::string JsonEscape(const std::string& value) {
  std::string output;
  output.reserve(value.size() + 8);
  for (const unsigned char c : value) {
    switch (c) {
      case '"': output += "\\\""; break;
      case '\\': output += "\\\\"; break;
      case '\b': output += "\\b"; break;
      case '\f': output += "\\f"; break;
      case '\n': output += "\\n"; break;
      case '\r': output += "\\r"; break;
      case '\t': output += "\\t"; break;
      default:
        if (c < 0x20) {
          const char hex[] = "0123456789abcdef";
          output += "\\u00";
          output += hex[c >> 4];
          output += hex[c & 0x0f];
        } else {
          output += static_cast<char>(c);
        }
    }
  }
  return output;
}

}  // namespace pgensparsescore
