cmake_minimum_required(VERSION 3.20)

project(pgensparsescore VERSION 0.1.0 LANGUAGES CXX)

include(FetchContent)
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PGENLIB_SOURCE_DIR "" CACHE PATH
    "Path to plink-ng/2.0/include; fetched from upstream when empty")

if(NOT PGENLIB_SOURCE_DIR)
  FetchContent_Declare(
    plink_ng
    GIT_REPOSITORY https://github.com/chrchang/plink-ng.git
    GIT_TAG 57b9e7d79ca51a7e66471927131af3a49ff5fa18
    GIT_SHALLOW FALSE
  )
  FetchContent_MakeAvailable(plink_ng)
  set(PGENLIB_SOURCE_DIR "${plink_ng_SOURCE_DIR}/2.0/include")
endif()

foreach(required_file
    pgenlib_read.cc
    pgenlib_ffi_support.cc
    plink2_base.cc
    plink2_bits.cc
    plink2_simd.cc
    pgenlib_misc.cc)
  if(NOT EXISTS "${PGENLIB_SOURCE_DIR}/${required_file}")
    message(FATAL_ERROR
      "PGENLIB_SOURCE_DIR does not look like plink-ng/2.0/include: "
      "missing ${required_file}")
  endif()
endforeach()

add_library(pgenlib STATIC
  "${PGENLIB_SOURCE_DIR}/plink2_base.cc"
  "${PGENLIB_SOURCE_DIR}/plink2_bits.cc"
  "${PGENLIB_SOURCE_DIR}/plink2_simd.cc"
  "${PGENLIB_SOURCE_DIR}/pgenlib_misc.cc"
  "${PGENLIB_SOURCE_DIR}/pgenlib_read.cc"
  "${PGENLIB_SOURCE_DIR}/pgenlib_ffi_support.cc"
)
target_include_directories(pgenlib PUBLIC "${PGENLIB_SOURCE_DIR}")
target_compile_options(pgenlib PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wno-unused-function>
)

find_package(ZLIB REQUIRED)

add_library(pgensparsescore_lib
  src/catalog.cc
  src/io.cc
  src/mapped_matrix.cc
  src/pgen_reader.cc
  src/scorer.cc
)
target_include_directories(pgensparsescore_lib PUBLIC src)
target_link_libraries(pgensparsescore_lib PUBLIC pgenlib ZLIB::ZLIB)
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
