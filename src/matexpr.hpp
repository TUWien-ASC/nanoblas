#ifndef FILE_MATEXPR
#define FILE_MATEXPR

#include <cstddef>
#include <iostream>

#include "vecexpr.hpp"

namespace nanoblas
{

  template <typename T>
  class MatExpr
  {
  public:
    auto upcast() const { return static_cast<const T&> (*this); }
    size_t rows() const { return upcast().rows(); }
    size_t cols() const { return upcast().cols(); }
    auto operator() (size_t i, size_t j) const { return upcast()(i,j); }
  };
  

  template <typename TM>
  std::ostream & operator<< (std::ostream & os, const MatExpr<TM> & m)
  {
    for (size_t i = 0; i < m.rows(); i++)
      {
        for (size_t j = 0; j < m.cols(); j++)
          os << m(i,j) << " ";
        os << std::endl;
      }
    return os;
  }

  template <typename TA, typename TB>
  class SumMatExpr : public MatExpr<SumMatExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    SumMatExpr (TA _a, TB _b) : a(_a), b(_b) { }
    auto operator() (size_t i) const { return a(i)+b(i); }
    size_t rows() const { return a.rows(); }
    size_t cols() const { return a.cols(); }      
  };
  
  template <typename TA, typename TB>
  auto operator+ (const MatExpr<TA> & a, const MatExpr<TB> & b)
  {
#ifndef NDEBUG
    if (a.rows() != b.rows() || a.cols() != b.cols())
      throw std::runtime_error("Matrix sizes do not match for addition, " 
                               "sizeof "+ std::to_string(a.rows())+"x"+std::to_string(a.cols()) + 
                               " and " + std::to_string(b.rows())+"x"+std::to_string(b.cols()) +    
                               " (in operator+).");     
#endif  
    return SumMatExpr(a.upcast(), b.upcast());
  }


  template <typename TA, typename TB>
  class MultMatMatExpr : public MatExpr<MultMatMatExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    MultMatMatExpr (TA _a, TB _b) : a(_a), b(_b) { }
    size_t rows() const { return a.rows(); }
    size_t cols() const { return b.cols(); }
    
    auto operator() (size_t i, size_t j) const { 
      // using elemtypeA = typename std::result_of<decltype(&operator())(TA,size_t,size_t)>::type;
      // using elemtypeB = typename std::result_of<decltype(&operator())(TB,size_t,size_t)>::type;
      using elemtypeA = typename std::invoke_result<decltype(&operator())(TA,size_t,size_t)>::type;
      using elemtypeB = typename std::invoke_result<decltype(&operator())(TB,size_t,size_t)>::type;
      using TSCAL = decltype(std::declval<elemtypeA>()*std::declval<elemtypeB>());
      
      TSCAL sum = 0;
      for (size_t k = 0; k < a.cols(); k++)
        sum += a(i,k) * b(k,j); 
      return sum;
    }
  };

  template <typename TA, typename TB>
  auto operator* (const MatExpr<TA> & a, const MatExpr<TB> & b)
  {
    return MultMatMatExpr<TA,TB>(a.upcast(), b.upcast());
  }
  
} // namespace nanoblas
#endif
