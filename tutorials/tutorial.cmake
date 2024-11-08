function(find_openmp FMMA_TARGET)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(${FMMA_TARGET} PRIVATE OpenMP::OpenMP_CXX)
  endif()
endfunction()

function(find_fmma TEST_TARGET)
  find_library(FMMA_LIBRARIES
    NAMES fmma
    /usr/local/lib
    REQUIRED
    )
  target_link_libraries(${TEST_TARGET} PRIVATE ${FMMA_LIBRARIES})
endfunction()

function(prepare TARGET_NAME test_source)
  add_executable(${TARGET_NAME} ${test_source})
  find_openmp(${TARGET_NAME})
  find_fmma(${TARGET_NAME})
  target_compile_options(${TARGET_NAME} PRIVATE -g -O3 -Wall -Wextra)
endfunction()
