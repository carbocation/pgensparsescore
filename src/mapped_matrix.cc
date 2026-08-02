// SPDX-License-Identifier: GPL-3.0-only
#include "mapped_matrix.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace pgensparsescore {

namespace {

std::runtime_error SystemError(const std::string& operation,
                               const std::string& path) {
  return std::runtime_error(operation + " " + path + ": " +
                            std::strerror(errno));
}

}  // namespace

MappedMatrix::MappedMatrix(const std::string& path, uint32_t row_ct,
                           uint32_t column_ct)
    : row_ct_(row_ct), column_ct_(column_ct) {
  if (!row_ct || !column_ct) {
    throw std::runtime_error("score matrix dimensions must be positive");
  }
  if (static_cast<uint64_t>(row_ct) >
      std::numeric_limits<uint64_t>::max() /
          (static_cast<uint64_t>(column_ct) * sizeof(double))) {
    throw std::runtime_error("score matrix size overflows uint64");
  }
  byte_ct_ = static_cast<uint64_t>(row_ct) * column_ct * sizeof(double);
  if (byte_ct_ > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
    throw std::runtime_error("score matrix is too large for this platform");
  }
  fd_ = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
  if (fd_ < 0) {
    throw SystemError("cannot create", path);
  }
  if (ftruncate(fd_, static_cast<off_t>(byte_ct_)) != 0) {
    const auto error = SystemError("cannot size", path);
    close(fd_);
    fd_ = -1;
    throw error;
  }
  void* mapping = mmap(nullptr, static_cast<size_t>(byte_ct_),
                       PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0);
  if (mapping == MAP_FAILED) {
    const auto error = SystemError("cannot memory-map", path);
    close(fd_);
    fd_ = -1;
    throw error;
  }
  data_ = static_cast<double*>(mapping);
}

MappedMatrix::~MappedMatrix() {
  if (data_) {
    msync(data_, static_cast<size_t>(byte_ct_), MS_SYNC);
    munmap(data_, static_cast<size_t>(byte_ct_));
  }
  if (fd_ >= 0) {
    close(fd_);
  }
}

double* MappedMatrix::Row(uint32_t row_idx) {
  if (row_idx >= row_ct_) {
    throw std::out_of_range("score matrix row out of range");
  }
  return data_ + static_cast<uint64_t>(row_idx) * column_ct_;
}

const double* MappedMatrix::Row(uint32_t row_idx) const {
  if (row_idx >= row_ct_) {
    throw std::out_of_range("score matrix row out of range");
  }
  return data_ + static_cast<uint64_t>(row_idx) * column_ct_;
}

void MappedMatrix::Flush() {
  if (msync(data_, static_cast<size_t>(byte_ct_), MS_SYNC) != 0) {
    throw std::runtime_error("cannot flush score matrix: " +
                             std::string(std::strerror(errno)));
  }
}

}  // namespace pgensparsescore
