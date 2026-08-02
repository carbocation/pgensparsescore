// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace pgensparsescore {

class MappedMatrix {
 public:
  MappedMatrix(const std::string& path, uint32_t row_ct, uint32_t column_ct);
  ~MappedMatrix();

  MappedMatrix(const MappedMatrix&) = delete;
  MappedMatrix& operator=(const MappedMatrix&) = delete;

  double* Row(uint32_t row_idx);
  const double* Row(uint32_t row_idx) const;
  uint32_t row_ct() const { return row_ct_; }
  uint32_t column_ct() const { return column_ct_; }
  uint64_t byte_ct() const { return byte_ct_; }
  void Flush();

 private:
  int fd_ = -1;
  double* data_ = nullptr;
  uint32_t row_ct_ = 0;
  uint32_t column_ct_ = 0;
  uint64_t byte_ct_ = 0;
};

}  // namespace pgensparsescore
