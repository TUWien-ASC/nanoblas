#include <iostream>

#include <matrix.hpp>
#include <lapack_interface.hpp>

using namespace nanoblas;

int main()
{
  Matrix<double> A(3, 3);
  Matrix<double> B(3, 3);
  Matrix<double> C(3, 3);

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      {
        A(i, j) = i + j;
        B(i, j) = i * j;
      }

  // C = A * B;
  MultMatMatLapack(A, B, C);
  
  std::cout << "A = \n" << A << std::endl;
  std::cout << "B = \n" << B << std::endl;
  std::cout << "C = A * B = \n" << C << std::endl;


  int n = 3;
  Matrix<double> a(n,n);
  for (int i = 0; i < n; i++)
    a(i,i) = 2;
  for (int i = 1; i < n; i++)
    a(i-1,i) = a(i,i-1) = -1;

  std::cout << "a = " << a << std::endl;

  std::cout << "Inv(a) = " << LapackLU(a).Inverse() << std::endl;
  
}
