cmake_minimum_required(VERSION 3.11)

project(fmma_test)

include(../cmake/common.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

function(set_test_properties TEST_TARGET)
  find_path(FMMA_INCLUDE_PATH
    NAMES fmma.hpp
    HINTS
    ${PROJECT_SOURCE_DIR}/../include/fmma
    REQUIRED
    )
  target_include_directories(${TEST_TARGET} PRIVATE ${FMMA_INCLUDE_PATH})

  find_library(FMMA_LIBRARIES
    NAMES fmma_cpu
    HINTS
    ${PROJECT_SOURCE_DIR}/../build
    REQUIRED
    )
  target_link_libraries(${TEST_TARGET} PRIVATE ${FMMA_LIBRARIES})
endfunction()

file(GLOB_RECURSE test_sources
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  src/*.cpp
  )

foreach(test_source ${test_sources})
  get_filename_component(MAIN_ITEM ${test_source} NAME_WE)
  get_filename_component(MAIN_DIRECTORY ${test_source} DIRECTORY)
  string(REPLACE "/" "_" COPY_DIRECTORY ${MAIN_DIRECTORY})
  string(CONCAT libname ${COPY_DIRECTORY} "_" ${MAIN_ITEM} "_cpu.out")
  set(TARGET_NAME ${libname})


  add_executable(${TARGET_NAME} ${test_source})
  set_common_properties(${TARGET_NAME})
  set_test_properties(${TARGET_NAME})
  target_compile_options(${TARGET_NAME} PRIVATE -g -O3 -Wall)
endforeach()
