# nanoblas core library

set(NANOBLAS_HEADERS
    vector.hpp
    vecexpr.hpp
    matrix.hpp
    matexpr.hpp
    lapack_interface.hpp
)


add_library(nanoblas INTERFACE)
target_include_directories(nanoblas INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(nanoblas INTERFACE cxx_std_20)

# Install headers
install(FILES ${NANOBLAS_HEADERS} DESTINATION nanoblas/include)

