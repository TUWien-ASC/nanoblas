cmake_minimum_required(VERSION 3.20)
project(nanoblas)

set (CMAKE_CXX_STANDARD 20)

include_directories(src)
include_directories(demos)

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NANOBIND_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NANOBIND_DIR}")
find_package(nanobind CONFIG REQUIRED)


if(WIN32)
   file (DOWNLOAD https://sourceforge.net/projects/openblas/files/v0.3.24/OpenBLAS-0.3.24-x64.zip/download "${CMAKE_BINARY_DIR}/openblas.zip" EXPECTED_MD5 fc08fe8c0dc7364da115d0e09b5a134f)
   file (ARCHIVE_EXTRACT INPUT ${CMAKE_BINARY_DIR}/openblas.zip DESTINATION ${CMAKE_BINARY_DIR}/openblas)
   set(BLAS_LIBRARIES ${CMAKE_BINARY_DIR}/openblas/lib/libopenblas.dll.a CACHE PATH "BLAS" FORCE)
   set(LAPACK_LIBRARIES ${CMAKE_BINARY_DIR}/openblas/lib/libopenblas.dll.a CACHE PATH "LAPACK" FORCE)
endif()


find_package(LAPACK REQUIRED)

add_subdirectory (demos)


nanobind_add_module(nanoblas_impl src/bind_bla.cpp)
target_link_libraries(nanoblas_impl PRIVATE LAPACK::LAPACK)

install (TARGETS nanoblas_impl DESTINATION nanoblas)
install (FILES src/vector.hpp DESTINATION nanoblas/include)


# Copy openblas dll file to build directory (no need to set PATH)
if(WIN32)
    add_custom_command(TARGET test_lapack POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${CMAKE_BINARY_DIR}/openblas/bin/libopenblas.dll"
      $<TARGET_FILE_DIR:test_lapack>
    )
endif()

