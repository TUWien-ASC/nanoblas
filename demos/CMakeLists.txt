add_executable (demo_vector demo_vector.cpp)
target_sources (demo_vector PUBLIC ../src/vector.hpp ../src/vecexpr.hpp)

add_executable (demo_matrix demo_matrix.cpp)
target_sources (demo_matrix PUBLIC ../src/matrix.hpp 
    ../src/vecexpr.hpp ../src/vector.hpp ../src/matexpr.hpp)
target_link_libraries (demo_matrix PUBLIC LAPACK::LAPACK)


add_executable (demo_lapack demo_lapack.cpp)
target_link_libraries (demo_lapack PUBLIC LAPACK::LAPACK)

