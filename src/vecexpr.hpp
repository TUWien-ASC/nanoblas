#ifndef FILE_EXPRESSION
#define FILE_EXPRESSION

namespace nanoblas
{

  template <typename T>
  class VecExpr
  {
  public:
    auto upcast() const { return static_cast<const T&> (*this); }
    size_t size() const { return upcast().size(); }
    auto operator() (size_t i) const { return upcast()(i); }
  };
  
 
  template <typename TA, typename TB>
  class SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    SumVecExpr (TA _a, TB _b) : a(_a), b(_b) { }

    auto operator() (size_t i) const { return a(i)+b(i); }
    size_t size() const { return a.size(); }      
  };
  
  template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA> & a, const VecExpr<TB> & b)
  {
#ifndef NDEBUG
    if (a.size() != b.size())
      throw std::runtime_error("Vector sizes do not match for addition" 
                               "sizeof "+ std::to_string(a.size()) + " and " + std::to_string(b.size()) +    
                               " (in operator+).");     
#endif  

    
    return SumVecExpr(a.upcast(), b.upcast());
  }



  
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
  
  template <typename T>
  auto operator* (double scal, const VecExpr<T> & v)
  {
    return ScaleVecExpr(scal, v.upcast());
  }



  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const VecExpr<T> & v)
  {
    if (v.size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.size(); i++)
      ost << ", " << v(i);
    return ost;
  }
  
}
 
#endif
