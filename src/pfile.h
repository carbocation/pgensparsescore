// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <string>
#include <vector>

namespace pgensparsescore {

struct PfileSpec {
  std::string pgen;
  std::string pvar;
  std::string psam;
};

std::vector<PfileSpec> ReadPfileList(const std::string& path);

}  // namespace pgensparsescore
