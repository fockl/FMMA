cmake_minimum_required(VERSION 3.11)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(fmma)

set(fmma_target "fmma_cpu")

include(cmake/common.cmake)
include(cmake/makelib.cmake)

