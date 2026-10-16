cmake_minimum_required(VERSION 3.18...3.31)
project(pydada2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(DADA2_ENABLE_COVERAGE "Build the C++ core with gcov instrumentation" OFF)
option(DADA2_BUILD_CPP_TESTS "Build the Catch2 C++ unit tests" OFF)
option(DADA2_FETCH_BOOST "Use the pinned standalone Boost.Math 1.86 (same on every platform); OFF = system headers" ON)

# SKBUILD_SABI_COMPONENT = Development.SABIModule when building abi3 wheels (wheel.py-api).
find_package(Python 3.11 COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT} REQUIRED)
find_package(nanobind CONFIG REQUIRED)

# ---------------------------------------------------------------- C++ core ---
file(GLOB_RECURSE DADA2_CORE_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/cpp/core/*.cpp")
add_library(dada2_core STATIC ${DADA2_CORE_SOURCES})
target_include_directories(dada2_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cpp/core")

# Floating-point reproducibility vs R: no FMA contraction, no fast-math.
if(MSVC)
  target_compile_options(dada2_core PRIVATE /fp:precise /W3)
  target_compile_definitions(dada2_core PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(dada2_core PRIVATE -ffp-contract=off -fno-fast-math -Wall -Wno-sign-compare
                         -Wno-unused-variable -Wno-unused-but-set-variable)
endif()

# Boost.Math (header-only) for ppois_upper (stats.cpp).
#   DADA2_FETCH_BOOST=ON (default): the standalone Boost.Math release pinned below, identical on
#     every platform (downloaded once at configure time; for offline builds pre-download it and
#     pass -DFETCHCONTENT_SOURCE_DIR_BOOST_MATH=<unpacked dir>).
#   DADA2_FETCH_BOOST=OFF: system/user headers only (never downloads). Point CMake at them with
#     -DDADA2_BOOST_MATH_INCLUDE=<dir containing boost/math> if they are not on the include path.
if(DADA2_FETCH_BOOST)
  message(STATUS "pydada2: using the pinned standalone Boost.Math 1.86")
  include(FetchContent)
  FetchContent_Declare(boost_math
    URL https://github.com/boostorg/math/archive/refs/tags/boost-1.86.0.tar.gz
    URL_HASH SHA256=c99287823d4b163e0c58df642557f1aaf2948e9b2f08211ff8c4f5ca1146ef8e
    SOURCE_SUBDIR do-not-add-subdirectory)  # headers only; don't configure its CMake project
  FetchContent_MakeAvailable(boost_math)
  target_include_directories(dada2_core SYSTEM PRIVATE "${boost_math_SOURCE_DIR}/include")
  target_compile_definitions(dada2_core PRIVATE BOOST_MATH_STANDALONE)
else()
  find_path(DADA2_BOOST_MATH_INCLUDE boost/math/special_functions/gamma.hpp)
  if(NOT DADA2_BOOST_MATH_INCLUDE)
    message(FATAL_ERROR
      "pydada2: DADA2_FETCH_BOOST=OFF but the Boost.Math headers "
      "(boost/math/special_functions/gamma.hpp) were not found. Install Boost (>= 1.86 recommended; "
      "older versions may change ppois_upper in the last bits), pass "
      "-DDADA2_BOOST_MATH_INCLUDE=<include dir>, or build with DADA2_FETCH_BOOST=ON.")
  endif()
  message(STATUS "pydada2: using Boost.Math headers from ${DADA2_BOOST_MATH_INCLUDE}")
  target_include_directories(dada2_core SYSTEM PRIVATE "${DADA2_BOOST_MATH_INCLUDE}")
  if(NOT EXISTS "${DADA2_BOOST_MATH_INCLUDE}/boost/config.hpp")  # standalone Boost.Math checkout
    target_compile_definitions(dada2_core PRIVATE BOOST_MATH_STANDALONE)
  endif()
endif()

find_package(Threads REQUIRED)
target_link_libraries(dada2_core PUBLIC Threads::Threads)

if(DADA2_ENABLE_COVERAGE AND NOT MSVC)
  target_compile_options(dada2_core PRIVATE --coverage -O0 -g)
  target_link_options(dada2_core PUBLIC --coverage)
endif()

# ----------------------------------------------------------- Python module ---
nanobind_add_module(_core NB_STATIC STABLE_ABI LTO
                    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/bindings/module.cpp")
target_link_libraries(_core PRIVATE dada2_core)
install(TARGETS _core LIBRARY DESTINATION pydada2)

# ----------------------------------------------------------- C++ unit tests ---
# cmake -DDADA2_BUILD_CPP_TESTS=ON ... && cmake --build <dir> && ctest --test-dir <dir>
if(DADA2_BUILD_CPP_TESTS)
  include(FetchContent)
  FetchContent_Declare(Catch2
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.7.1.tar.gz
    URL_HASH SHA256=c991b247a1a0d7bb9c39aa35faf0fe9e19764213f28ffba3109388e62ee0269c)
  FetchContent_MakeAvailable(Catch2)
  file(GLOB DADA2_CPP_TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests/*.cpp")
  add_executable(dada2_cpp_tests ${DADA2_CPP_TESTS})
  target_link_libraries(dada2_cpp_tests PRIVATE dada2_core Catch2::Catch2WithMain)
  enable_testing()
  add_test(NAME dada2_cpp_tests COMMAND dada2_cpp_tests)
endif()
