
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

set_target_properties(${fmma_target} PROPERTIES OUTPUT_NAME ${fmma_link})

install(
  TARGETS ${fmma_target}
  LIBRARY
  DESTINATION lib
  )

