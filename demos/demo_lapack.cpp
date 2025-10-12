#include <iostream>

#include <vector.hpp>
#include <lapack_interface.hpp>


using namespace nanoblas;
using namespace std;


int main()
{
  Vector<double> x(5);
  Vector<double> y(5);

  for (int i = 0; i < x.size(); i++)
    {
      x(i) = i;
      y(i) = 2;
    }

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  
  AddVectorLapack (2, x, y);  
  cout << "y+2*x = " << y << endl;

  Matrix<double> A(5,5);
  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      A(i,j) = i + j + 1; 

  MultMatVecLapack<double> (1.0, A, x, 0.0, y);
  cout << "A*x = " << y << endl;
  cout << " =?=" << A * x << endl;
}
