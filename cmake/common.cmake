function(set_common_properties FMMA_TARGET)
  target_compile_options(${FMMA_TARGET} PRIVATE -g -O3 -Wall -Wextra -fPIC)

  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    add_definitions(-DOPENMP_FOUND)
    message(STATUS "Found OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    target_compile_options(${FMMA_TARGET} PRIVATE -fopenmp)
    target_link_libraries(${FMMA_TARGET} PRIVATE OpenMP::OpenMP_CXX)
  else()
    message(STATUS "Not Found OpenMP")
  endif()

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
