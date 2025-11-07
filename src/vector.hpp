#ifndef FILE_VECTOR
#define FILE_VECTOR

#include <iostream>
#include <vector>

#include "vecexpr.hpp"


namespace nanoblas
{
 
  enum ORDERING { RowMajor, ColMajor };

  template <typename T, ORDERING ORD = RowMajor>
  class MatrixView;


  
  template <typename T=double, typename TDIST = std::integral_constant<size_t,1> >
  class VectorView : public VecExpr<VectorView<T,TDIST>>
  {
  protected:
    T* m_data;
    size_t m_size;
    TDIST m_dist;
    
  public:
    VectorView() = default;
    VectorView(const VectorView&) = default;
    
    template <typename TDIST2>
    VectorView (const VectorView<T,TDIST2>& v2)
      : m_data(v2.data()), m_size(v2.size()), m_dist(v2.dist()) { }
    
    VectorView (size_t size, T* data)
      : m_data(data), m_size(size) { }
    
    VectorView (size_t size, TDIST dist, T* data)
      : m_data(data), m_size(size), m_dist(dist) { }
    
    VectorView operator= (const VectorView& v2)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = v2(i);
      return *this;
    }

    template <typename TB>
    VectorView operator= (const VecExpr<TB>& v2)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = v2(i);
      return *this;
    }

    VectorView operator= (T scal)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = scal;
      return *this;
    }

    VectorView operator= (const std::vector<T>& v2)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = v2[i];
      return *this;
    }

    operator std::vector<T>() const
    {
      std::vector<T> v2(m_size);
      for (size_t i = 0; i < m_size; i++)
        v2[i] = m_data[m_dist*i];
      return v2;
    }

    T * data() const { return m_data; }
    size_t size() const { return m_size; }
    auto dist() const { return m_dist; }
    
    T& operator()(size_t i) { return m_data[m_dist*i]; }
    const T& operator()(size_t i) const { return m_data[m_dist*i]; }
    
    T& operator[](size_t i) { return m_data[m_dist*i]; }
    const T& operator[](size_t i) const { return m_data[m_dist*i]; }
    
    auto range(size_t first, size_t next) const {
      return VectorView(next-first, m_dist, m_data+first*m_dist);
    }

    auto slice(size_t first, size_t slice) const {
      return VectorView<T,size_t> (m_size/slice, m_dist*slice, m_data+first*m_dist);
    }

    auto asMatrix(size_t rows, size_t cols) const {
      return MatrixView<T,RowMajor> (rows, cols, m_data);  
    }
    
    template <typename TB>
    VectorView& operator+= (const VecExpr<TB>& v2)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] += v2(i);
      return *this;
    }

    template <typename TB>
    VectorView& operator-= (const VecExpr<TB>& v2)
      {
        for (size_t i = 0; i < m_size; i++)
          m_data[m_dist*i] -= v2(i);
        return *this;
      }

    VectorView& operator*= (T scal)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] *= scal;
      return *this;
    }
    
  };
  
  

  
  template <typename T=double>
  class Vector : public VectorView<T>
  {
    typedef VectorView<T> BASE;
    using BASE::m_size;
    using BASE::m_data;
  public:
    explicit Vector (size_t size) 
      : VectorView<T> (size, new T[size]) { ; }
    
    Vector (const Vector& v)
      : Vector(v.size())
    {
      *this = v;
    }

    Vector (Vector && v)
      : VectorView<T> (0, nullptr)
    {
      std::swap(m_size, v.m_size);
      std::swap(m_data, v.m_data);
    }

    template <typename TB>
    Vector (const VecExpr<TB>& v)
      : Vector(v.size())
    {
      *this = v;
    }

  
    Vector (std::initializer_list<T> list) 
      : VectorView<T> (list.size(), new T[list.size()])
    {
      size_t cnt = 0;
      for (auto val : list)
        (*this)(cnt++) = val;
    }
    
    ~Vector () { delete [] m_data; }

    using BASE::operator=;
    Vector& operator=(const Vector& v2)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[i] = v2(i);
      return *this;
    }

    Vector& operator= (Vector && v2)
    {
      std::swap(m_size, v2.m_size);
      std::swap(m_data, v2.m_data);
      return *this;
    }
  };


  template <typename ...Args>
  std::ostream& operator<< (std::ostream& ost, const VectorView<Args...>& v)
  {
    if (v.size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.size(); i++)
      ost << ", " << v(i);
    return ost;
  }


  template <size_t S, typename T=double>
  class Vec : public VecExpr<Vec<S,T>>
  {
    std::array<T, S> m_data;
  public:
    Vec() = default;
    Vec (const Vec& v) = default;
    
    template <typename TB>
    Vec (const VecExpr<TB>& v)
    {
      for (size_t i = 0; i < S; i++)
        m_data[i] = v(i);
    }

    Vec (T val)
    {
      for (size_t i = 0; i < S; i++)
        m_data[i] = val;
    }

    Vec (std::initializer_list<T> list)
    {
      size_t cnt = 0;
      for (auto val : list)
        (*this)(cnt++) = val;
    }
    
    Vec& operator= (const Vec& v2)
    {
      for (size_t i = 0; i < S; i++)
        m_data[i] = v2.m_data[i];
      return *this;
    }
    
    template <typename TB>
    Vec& operator= (const VecExpr<TB>& v2)
    {
      for (size_t i = 0; i < S; i++)
        m_data[i] = v2(i);
      return *this;
    }
    
    size_t size() const { return S; }
    std::array<T,S>& data() { return m_data; }
    const std::array<T,S>& data() const { return m_data; }
    
    T& operator() (size_t i) { return m_data[i]; }
    const T& operator() (size_t i) const { return m_data[i]; }
 };
  
}

#endif
