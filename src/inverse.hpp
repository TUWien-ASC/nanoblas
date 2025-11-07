#ifndef INVERSE_HPP
#define INVERSE_HPP

#include <vector>

#include <matrix.hpp>

namespace nanoblas {

template <typename T>  
void calcInverse(const Matrix<T>& mat) 
{
    // Check if the matrix is square
    if (mat.rows() != mat.cols()) {
        throw std::invalid_argument("Matrix must be square to compute its inverse.");
    }

    size_t n = mat.rows();

    std::vector<int> p(n);   // pivot-permutation
    for (size_t j = 0; j < n; j++) p[j] = j;

    for (size_t j = 0; j < n; j++)
      {
	// pivot search
	double maxval = abs(mat(j,j));
	size_t r = j;

	for (size_t i = j+1; i < n; i++)
	  if (abs (mat(j, i)) > maxval)
	    {
	      r = i;
	      maxval = abs (mat(r, i));
	    }
      
        double rest = 0.0;
        for (size_t i = j+1; i < n; i++)
          rest += abs(mat(r, i));
	if (maxval < 1e-20*rest)
	  {
	    throw std::runtime_error("Inverse matrix: Matrix singular");
	  }

	// exchange rows
	if (r > j)
	  {
	    for (size_t k = 0; k < n; k++)
	      std::swap (mat(k, j), mat(k, r));
	    std::swap (p[j], p[r]);
	  }
      

	// transformation
	
	T hr = 1.0 / mat(j,j);

	for (size_t i = 0; i < n; i++)
            mat(j,i) = hr * mat(j,i);
	mat(j,j) = hr;

	for (size_t k = 0; k < n; k++)
	  if (k != j)
	    {
	      T help = mat(n*k+j);
	      T h = help * hr;   

	      for (size_t i = 0; i < n; i++)
		{
		  T h = help * mat(n*j+i); 
		  mat(n*k+i) -= h;
		}

	      mat(k,j) = -h;
	    }
      }

    // row exchange
  
    std::vector<T> hv(n);
    for (size_t i = 0; i < n; i++)
      {
	for (size_t k = 0; k < n; k++) hv[p[k]] = mat(k, i);
	for (size_t k = 0; k < n; k++) mat(k, i) = hv[k];
      }

    }


}

#endif


