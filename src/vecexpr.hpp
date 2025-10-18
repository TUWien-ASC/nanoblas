#ifndef FILE_EXPRESSION
#define FILE_EXPRESSION

#include<complex>
#include<cassert>
#include <type_traits>

namespace nanoblas
{

  

  /*
    Expression templates for vector expressions
    
    VecExpr is the basis for all vector expressions
    Uses compile-time polymorphism, known as
    Curiously recurring template pattern (CRTP)
  */

  template <typename T>
  class VecExpr
  {
  public:
    auto derived() const { return static_cast<const T&> (*this); }
    size_t size() const { return derived().size(); }
    auto operator() (size_t i) const { return derived()(i); }
  };

  

  // ************************ SumVecExpr *********************
 
  template <typename TA, typename TB>
  struct SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    SumVecExpr (TA _a, TB _b) : a(_a), b(_b) { }

    auto operator() (size_t i) const { return a(i)+b(i); }
    size_t size() const { return a.size(); }      
  };
  
  template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA>& a, const VecExpr<TB>& b)
  {
    assert(a.size()==b.size());
    return SumVecExpr(a.derived(), b.derived());
  }


  // ************************ SubVecExpr *********************
  
  template <typename TA, typename TB>
  class SubVecExpr : public VecExpr<SubVecExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    SubVecExpr (TA _a, TB _b) : a(_a), b(_b) { }

    auto operator() (size_t i) const { return a(i)-b(i); }
    size_t size() const { return a.size(); }      
  };
  
  template <typename TA, typename TB>
  auto operator- (const VecExpr<TA>& a, const VecExpr<TB>& b)
  {
    assert(a.size()==b.size());    
    return SubVecExpr(a.derived(), b.derived());
  }


  // ************************ NegVecExpr *********************

  
  template <typename TA>
  class NegVecExpr : public VecExpr<NegVecExpr<TA>>
  {
    TA a;
  public:
    NegVecExpr (TA _a) : a(_a) { }

    auto operator() (size_t i) const { return -a(i); }
    size_t size() const { return a.size(); }      
  };
  
  template <typename TA>
  auto operator- (const VecExpr<TA>& a)
  {
    return NegVecExpr(a.derived());
  }



  // ************************ ScaleVecExpr *********************
  

  template <typename T>
  struct is_scalar_type { static constexpr bool value = std::integral<T>||std::floating_point<T>; };
  
  template <typename T>
  struct is_scalar_type<std::complex<T>> { static constexpr bool value = isScalar<T>(); };


  template <typename T>
  constexpr bool isScalar() { return is_scalar_type<T>::value; }


  
  template <typename TSCAL, typename TV>
  class ScaleVecExpr : public VecExpr<ScaleVecExpr<TSCAL,TV>>
  {
    TSCAL scal;
    TV vec;
  public:
    ScaleVecExpr (TSCAL _scal, TV _vec) : scal(_scal), vec(_vec) { }
    auto operator() (size_t i) const { return scal*vec(i); }
    size_t size() const { return vec.size(); }      
  };

  /*
  // C++17 needs the enable_if workaround instead of the 'requires' concept syntax
  template <typename TSCAL, typename T,
            typename std::enable_if<isScalar<TSCAL>(),int>::type = 0>
  */
  
  template <typename TSCAL, typename T> requires (isScalar<TSCAL>())
  auto operator* (TSCAL scal, const VecExpr<T>& v)
  {
    return ScaleVecExpr(scal, v.derived());
  }


  // **************** dot product of two vectors *****************
 
  template <typename TA, typename TB>
  auto dot (const VecExpr<TA>& a, const VecExpr<TB>& b)
  {
    assert (a.size() == b.size());

    using elemtypeA = typename std::invoke_result<TA,size_t>::type;
    using elemtypeB = typename std::invoke_result<TB,size_t>::type;
    using TSUM = decltype(std::declval<elemtypeA>()*std::declval<elemtypeB>());

    TSUM sum = 0;
    for (size_t i = 0; i < a.size(); i++)
      sum += a(i)*b(i);
    return sum;
  }


  
  // **************** euclidean norm of vector *****************

  double norm2 (double x) { return x*x; }
  double norm2 (std::complex<double> x) { return x.real()*x.real() + x.imag()*x.imag(); }

  
  template <typename TA>
  auto norm (const VecExpr<TA>& a)
  {
    using elemtype = typename std::remove_cvref<typename std::invoke_result<TA,size_t>::type>::type;

    elemtype sum = 0;
    for (size_t i = 0; i < a.size(); i++)
      sum += norm2(a(i));
    return sqrt(sum);
  }

  
  
  // ***********************  output operator  *********************
  

  template <typename T>
  std::ostream & operator<< (std::ostream& ost, const VecExpr<T>& v)
  {
    if (v.size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.size(); i++)
      ost << ", " << v(i);
    return ost;
  }
  
}
 
#endif
