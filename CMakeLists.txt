cmake_minimum_required(VERSION 3.11)

set(CMAKE_OPTIMIZE_DEPENDENCIES ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(BUILD_SHARED_LIBS ON)

project(fmma)

set(fmma_target "fmma_cpu")

option(FMMA_USE_BLAS "Use BLAS" OFF)

include(cmake/common.cmake)
include(cmake/makelib.cmake)

