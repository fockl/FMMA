function(find_openmp FMMA_TARGET)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(${FMMA_TARGET} PRIVATE OpenMP::OpenMP_CXX)
  endif()
endfunction()

function(find_fmma TEST_TARGET)
  find_library(FMMA_LIBRARIES
    NAMES fmma
    HINTS
    /usr/local/lib
    ${CMAKE_CURRENT_LIST_DIR}/../../build
    REQUIRED
    )
  target_link_libraries(${TEST_TARGET} PRIVATE ${FMMA_LIBRARIES})

  find_path(FMMA_INCLUDE_PATH
    NAMES fmma
    HINTS
    /usr/local/include
    ${CMAKE_CURRENT_LIST_DIR}/../../include
    REQUIRED
    )
  target_include_directories(${TEST_TARGET} PRIVATE ${FMMA_INCLUDE_PATH})
endfunction()

function(prepare TARGET_NAME test_source)
  add_executable(${TARGET_NAME} ${test_source})
  find_openmp(${TARGET_NAME})
  find_fmma(${TARGET_NAME})
  target_compile_options(${TARGET_NAME} PRIVATE -g -O3 -Wall -Wextra)
endfunction()
