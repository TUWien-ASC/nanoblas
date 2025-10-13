#ifndef FILE_MATEXPR
#define FILE_MATEXPR

#include <cstddef>
#include <iostream>

#include "vecexpr.hpp"

namespace nanoblas
{

  /*
    Expression templates for matrix expressions
    and matrix-vector expressions
  */
  
  
  // base class for all matrix expressions
  
  template <typename T>
  class MatExpr
  {
  public:
    auto upcast() const { return static_cast<const T&> (*this); }
    size_t rows() const { return upcast().rows(); }
    size_t cols() const { return upcast().cols(); }
    auto operator() (size_t i, size_t j) const { return upcast()(i,j); }
  };
  
  // ************************* output operator *******************
  
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


  // ************************* SumMatExpr *******************  

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
    assert(a.rows()==b.rows() && a.cols()==b.cols());
    return SumMatExpr(a.upcast(), b.upcast());
  }

  

  // ************************* MultMatMatExpr *******************
  
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
      using elemtypeA = std::invoke_result<TA,size_t,size_t>::type;
      using elemtypeB = std::invoke_result<TB,size_t,size_t>::type;
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
    assert(a.cols()==b.rows());
    return MultMatMatExpr<TA,TB>(a.upcast(), b.upcast());
  }
  


  // ************************* MultMatVecExpr *******************

 template <typename TA, typename TB>
  class MultMatVecExpr : public VecExpr<MultMatVecExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    MultMatVecExpr (TA _a, TB _b) : a(_a), b(_b) { }
    size_t size() const { return a.rows(); }
    
    auto operator() (size_t i) const { 
      using elemtypeA = std::invoke_result<TA,size_t,size_t>::type;
      using elemtypeB = std::invoke_result<TB,size_t>::type;
      using TSCAL = decltype(std::declval<elemtypeA>()*std::declval<elemtypeB>());
      
      TSCAL sum = 0;
      for (size_t k = 0; k < a.cols(); k++)
        sum += a(i,k) * b(k); 
      return sum;
    }
  };

  template <typename TA, typename TB>
  auto operator* (const MatExpr<TA> & a, const VecExpr<TB> & b)
  {
    assert(a.cols()==b.size());    
    return MultMatVecExpr<TA,TB>(a.upcast(), b.upcast());
  }
 



} // namespace nanoblas

#endif
