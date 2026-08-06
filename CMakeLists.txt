cmake_minimum_required(VERSION 3.20)

project(pgensparsescore VERSION 0.1.0 LANGUAGES CXX)

include(FetchContent)
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PLINK_NG_SOURCE_DIR "" CACHE PATH
    "Path to a plink-ng checkout; the pinned repository is fetched when empty")
set(PGENLIB_SOURCE_DIR "" CACHE PATH
    "Legacy path to plink-ng/2.0/include")

if(PLINK_NG_SOURCE_DIR AND PGENLIB_SOURCE_DIR)
  message(FATAL_ERROR
    "Set only one of PLINK_NG_SOURCE_DIR and PGENLIB_SOURCE_DIR")
endif()

set(PGEN_RANS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(PGEN_RANS_INSTALL OFF CACHE BOOL "" FORCE)
if(PLINK_NG_SOURCE_DIR OR PGENLIB_SOURCE_DIR)
  if(PGENLIB_SOURCE_DIR)
    get_filename_component(_plink_ng_2_dir "${PGENLIB_SOURCE_DIR}" DIRECTORY)
    get_filename_component(_plink_ng_source_dir "${_plink_ng_2_dir}" DIRECTORY)
  else()
    set(_plink_ng_source_dir "${PLINK_NG_SOURCE_DIR}")
  endif()
  if(NOT EXISTS "${_plink_ng_source_dir}/2.0/pgen_rans_pgenlib.h")
    message(FATAL_ERROR
      "The configured plink-ng checkout does not export the rANS pgenlib adapter")
  endif()
  add_subdirectory(
    "${_plink_ng_source_dir}"
    "${CMAKE_CURRENT_BINARY_DIR}/plink_ng"
  )
else()
  FetchContent_Declare(
    plink_ng
    GIT_REPOSITORY https://github.com/carbocation/plink-ng.git
    GIT_TAG 9b3dc0c4397e239214c441ad170a567b1138313c
    GIT_SHALLOW FALSE
  )
  FetchContent_MakeAvailable(plink_ng)
endif()

find_package(ZLIB REQUIRED)
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h REQUIRED)
find_library(ZSTD_LIBRARY NAMES zstd libzstd REQUIRED)
add_library(ZSTD::ZSTD UNKNOWN IMPORTED)
set_target_properties(ZSTD::ZSTD PROPERTIES
  IMPORTED_LOCATION "${ZSTD_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${ZSTD_INCLUDE_DIR}"
)

add_library(pgensparsescore_lib
  src/catalog.cc
  src/compiled_catalog.cc
  src/frequency.cc
  src/fragment_scorer.cc
  src/io.cc
  src/mapped_matrix.cc
  src/pfile.cc
  src/pgen_reader.cc
  src/progress.cc
  src/scorer.cc
  src/score_fragment.cc
  src/variant_index.cc
)
target_include_directories(pgensparsescore_lib PUBLIC src)
target_link_libraries(
  pgensparsescore_lib PUBLIC pgen_rans::pgenlib ZLIB::ZLIB ZSTD::ZSTD
)
target_compile_options(pgensparsescore_lib PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall;-Wextra;-Wpedantic>
)

add_executable(pgensparsescore src/main.cc)
target_link_libraries(pgensparsescore PRIVATE pgensparsescore_lib)

if(BUILD_TESTING)
  add_executable(unit_tests tests/unit_tests.cc)
  target_link_libraries(unit_tests PRIVATE pgensparsescore_lib)
  add_test(NAME unit_tests COMMAND unit_tests)

  find_package(Python3 COMPONENTS Interpreter)
  find_program(PLINK2_EXECUTABLE NAMES plink2)
  if(Python3_Interpreter_FOUND)
    add_test(
      NAME progress_integration
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/tests/progress_test.py"
              --scorer "$<TARGET_FILE:pgensparsescore>"
    )
  endif()
  if(Python3_Interpreter_FOUND AND PLINK2_EXECUTABLE)
    add_test(
      NAME pgen_integration
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration_test.py"
              --plink2 "${PLINK2_EXECUTABLE}"
              --scorer "$<TARGET_FILE:pgensparsescore>"
    )
  else()
    message(STATUS "plink2 not found; skipping real-PGEN integration test")
  endif()
endif()

install(TARGETS pgensparsescore RUNTIME DESTINATION bin)
