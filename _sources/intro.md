# Welcome to nanoblas' documentation!


nanoblas is a C++ library for basic linear algebra operations.
The library provides template classes **Vector** and **Matrix**.

## Installation

install it via git-clone:

    git clone https://github.com/TUWien-ASC/nanoblas.git


To configure and build some tests do

    cd nanoblas
    mkdir build
    cd build
    cmake ..
    make
    

## Using nanoblas

To use ASC-bla in your code, set the compiler include path properly, and include the header files

    #include <vector.hpp>
    #include <matrix.hpp>

All objects are implemented in the namespace ASC_bla. To use them with less typing, you can set

    namespace bla = nanoblas;

or even

    
    using namespace nanoblas;

    

You can create vectors and compute with vectors like:

                 
```cpp
Vector<double> x(5), y(5), z(5);
for (int i = 0; i < x.size(); i++)
   x(i) = i;
y = 5.0
z = x+3*y;
cout << "z = " << z << endl;
```

For matrices you can choose between row-major (`RowMajor`) or column-major (`ColMajor`) storage,
default is row-major.

```cpp
Matrix<double,RowMajor> m1(5,3), m2(3,3);
for (int i = 0; i < m1.rows(); i++)
  for (int j = 0; j < m1.cols(); j++)
    m1(i,j) = i+j;
m2 = 3.7;
Matrix<double> product = m1 * m2;
```

You can extract a row or a column of a matrix:

```cpp
Vector col1 = product.col(1);
```

some changes ...  

   
