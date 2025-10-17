# Demos for nanoblas

set(NANOBLAS_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src")

# Demo: demo_vector
add_executable(demo_vector demo_vector.cpp)
target_include_directories(demo_vector PRIVATE "${NANOBLAS_SRC_DIR}")
target_compile_features(demo_vector PRIVATE cxx_std_20)

# Demo: demo_matrix
add_executable(demo_matrix demo_matrix.cpp)
target_include_directories(demo_matrix PRIVATE "${NANOBLAS_SRC_DIR}")
target_link_libraries(demo_matrix PRIVATE LAPACK::LAPACK)
target_compile_features(demo_matrix PRIVATE cxx_std_20)

# Demo: demo_lapack
add_executable(demo_lapack demo_lapack.cpp)
target_include_directories(demo_lapack PRIVATE "${NANOBLAS_SRC_DIR}")
target_link_libraries(demo_lapack PRIVATE LAPACK::LAPACK)
target_compile_features(demo_lapack PRIVATE cxx_std_20)

# Install demo executables (optional)
install(TARGETS demo_vector demo_matrix demo_lapack
    RUNTIME DESTINATION nanoblas/demo
)