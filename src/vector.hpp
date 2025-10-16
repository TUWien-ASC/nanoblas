#ifndef FILE_VECTOR
#define FILE_VECTOR

#include <iostream>

#include "vecexpr.hpp"


namespace nanoblas
{
 
  template <typename T, typename TDIST = std::integral_constant<size_t,1> >
  class VectorView : public VecExpr<VectorView<T,TDIST>>
  {
  protected:
    T * data_;
    size_t size_;
    TDIST dist_;
  public:
    VectorView() = default;
    VectorView(const VectorView &) = default;
    
    template <typename TDIST2>
    VectorView (const VectorView<T,TDIST2> & v2)
      : data_(v2.data()), size_(v2.size()), dist_(v2.dist()) { }
    
    VectorView (size_t size, T * data)
      : data_(data), size_(size) { }
    
    VectorView (size_t size, TDIST dist, T * data)
      : data_(data), size_(size), dist_(dist) { }
    
    template <typename TB>
    VectorView & operator= (const VecExpr<TB> & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] = v2(i);
      return *this;
    }

    VectorView & operator= (T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] = scal;
      return *this;
    }

    T * data() const { return data_; }
    size_t size() const { return size_; }
    auto dist() const { return dist_; }
    
    T & operator()(size_t i) { return data_[dist_*i]; }
    const T & operator()(size_t i) const { return data_[dist_*i]; }
    
    auto range(size_t first, size_t next) const {
      return VectorView(next-first, dist_, data_+first*dist_);
    }

    auto slice(size_t first, size_t slice) const {
      return VectorView<T,size_t> (size_/slice, dist_*slice, data_+first*dist_);
    }

    template <typename TB>
    VectorView & operator+= (const VecExpr<TB> & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] += v2(i);
      return *this;
    }

    template <typename TB>
    VectorView & operator-= (const VecExpr<TB> & v2)
      {
        for (size_t i = 0; i < size_; i++)
          data_[dist_*i] -= v2(i);
        return *this;
      }
  };
  
  

  
  template <typename T=double>
  class Vector : public VectorView<T>
  {
    typedef VectorView<T> BASE;
    using BASE::size_;
    using BASE::data_;
  public:
    Vector (size_t size) 
      : VectorView<T> (size, new T[size]) { ; }
    
    Vector (const Vector & v)
      : Vector(v.size())
    {
      *this = v;
    }

    Vector (Vector && v)
      : VectorView<T> (0, nullptr)
    {
      std::swap(size_, v.size_);
      std::swap(data_, v.data_);
    }

    template <typename TB>
    Vector (const VecExpr<TB> & v)
      : Vector(v.size())
    {
      *this = v;
    }
    
    ~Vector () { delete [] data_; }

    using BASE::operator=;
    Vector & operator=(const Vector & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = v2(i);
      return *this;
    }

    Vector & operator= (Vector && v2)
    {
      std::swap(size_, v2.size_);
      std::swap(data_, v2.data_);
      return *this;
    }
  };


  template <typename ...Args>
  std::ostream & operator<< (std::ostream & ost, const VectorView<Args...> & v)
  {
    if (v.size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.size(); i++)
      ost << ", " << v(i);
    return ost;
  }
  
}

#endif
