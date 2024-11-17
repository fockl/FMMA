
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE fmma_sources
  RELATIVE ${PROJECT_SOURCE_DIR}
  src/*.cpp
  )

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
  DESTINATION include
  FILES_MATCHING PATTERN *.hpp
  )

add_library(${fmma_target} ${fmma_sources})
set_common_properties(${fmma_target})
set_blas(${fmma_target})

if(FMMA_PYTHON)
  add_subdirectory(extern/pybind11)
  pybind11_add_module(pyfmma cpython/pyfmma.cpp)
  target_link_libraries(pyfmma PRIVATE ${fmma_target})
endif()

set_target_properties(${fmma_target} PROPERTIES OUTPUT_NAME ${fmma_link})

install(
  TARGETS ${fmma_target}
  LIBRARY
  DESTINATION lib
  )

