cmake_minimum_required(VERSION 3.20)
project(nanoblas LANGUAGES CXX)

option(NANOBLAS_BUILD_DEMOS "Build demonstration targets" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Dependency: nanobind via Python
find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NANOBIND_DIR
)
list(APPEND CMAKE_PREFIX_PATH "${NANOBIND_DIR}")
find_package(nanobind CONFIG REQUIRED)

# BLAS/LAPACK handling (Windows-specific OpenBLAS handling)
if(WIN32)
    set(OPENBLAS_URL "https://sourceforge.net/projects/openblas/files/v0.3.24/OpenBLAS-0.3.24-x64.zip/download")
    set(OPENBLAS_ZIP "${CMAKE_BINARY_DIR}/openblas.zip")
    set(OPENBLAS_MD5 "fc08fe8c0dc7364da115d0e09b5a134f")

    file(DOWNLOAD ${OPENBLAS_URL} ${OPENBLAS_ZIP} EXPECTED_MD5 ${OPENBLAS_MD5})
    file(ARCHIVE_EXTRACT INPUT ${OPENBLAS_ZIP} DESTINATION ${CMAKE_BINARY_DIR}/openblas)

    set(BLAS_LIBRARIES "${CMAKE_BINARY_DIR}/openblas/lib/libopenblas.dll.a" CACHE PATH "BLAS" FORCE)
    set(LAPACK_LIBRARIES "${CMAKE_BINARY_DIR}/openblas/lib/libopenblas.dll.a" CACHE PATH "LAPACK" FORCE)
endif()

find_package(LAPACK REQUIRED)

# Create an interface library that carries the LAPACK dependency
add_library(nanoblas INTERFACE)
target_link_libraries(nanoblas INTERFACE ${LAPACK_LIBRARIES})
target_include_directories(nanoblas INTERFACE ${LAPACK_INCLUDE_DIRS})



nanobind_add_module(nanoblas_impl src/bind_bla.cpp)
target_include_directories(nanoblas_impl PRIVATE src)
target_link_libraries(nanoblas_impl PRIVATE LAPACK::LAPACK)
target_compile_features(nanoblas_impl PRIVATE cxx_std_20)

install(TARGETS nanoblas_impl DESTINATION nanoblas)
install(FILES src/vector.hpp DESTINATION nanoblas/include)

# Demos (conditional)
if(NANOBLAS_BUILD_DEMOS)
    add_subdirectory(demos)
endif()

# Windows: copy openblas DLL for demo_lapack after build
if(WIN32 AND TARGET demo_lapack)
    add_custom_command(TARGET demo_lapack POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_BINARY_DIR}/openblas/bin/libopenblas.dll"
            $<TARGET_FILE_DIR:demo_lapack>
    )
endif()
