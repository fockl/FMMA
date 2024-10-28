function(set_blas FMMA_TARGET)
  if(FMMA_USE_BLAS)
    target_compile_definitions(${FMMA_TARGET} PRIVATE FMMA_USE_BLAS)
    message(STATUS "Use Blas")
    find_package(BLAS REQUIRED)
    find_path(BLAS_INCLUDE_DIRS
      NAMES cblas.h
      HINTS
      /usr/include
      /usr/local/include
      /usr/include/openblas
      )
    message(STATUS ${BLAS_INCLUDE_DIRS})

    find_library(BLAS_LIBRARIES_NEW
      NAMES blas
      HINTS
      /usr/lib
      /usr/local/lib
      )

    target_include_directories(${FMMA_TARGET} PRIVATE ${BLAS_INCLUDE_DIRS})
    target_link_libraries(${FMMA_TARGET} PRIVATE ${BLAS_LIBRARIES} ${BLAS_LIBRARIES})
    target_link_libraries(${FMMA_TARGET} PRIVATE ${BLAS_LIBRARIES} ${BLAS_LIBRARIES_NEW})

    message(STATUS ${BLAS_LIBRARIES})
    message(STATUS ${BLAS_LIBRARIES_NEW})
  else()
    message(STATUS "Not Use Blas")
  endif()

endfunction()
