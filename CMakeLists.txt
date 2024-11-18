cmake_minimum_required(VERSION 3.12)

set(CMAKE_OPTIMIZE_DEPENDENCIES ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(BUILD_SHARED_LIBS ON)

project(fmma)

set(fmma_target "fmma")
set(fmma_link "fmma")

option(FMMA_CMAKE_DEBUG_FLAG "Cmake Debug flag" ON)
option(FMMA_USE_BLAS "Use BLAS" OFF)
option(FMMA_DEBUG_MODE "Debug mode" OFF)
option(FMMA_TIME_LOG "Time log" OFF)
option(FMMA_PYTHON "Create Python library" OFF)

include(cmake/set_name.cmake)
include(cmake/options.cmake)
include(cmake/common.cmake)
include(cmake/blas.cmake)
include(cmake/gpu.cmake)
include(cmake/makelib.cmake)

