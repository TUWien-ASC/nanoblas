#ifndef FILE_MATRIX
#define FILE_MATRIX

#include "matexpr.hpp"
#include "vector.hpp"

namespace nanoblas
{
  
  // enum ORDERING { RowMajor, ColMajor };

  template <typename T, ORDERING ORD>
  class MatrixView : public MatExpr<MatrixView<T,ORD>>
  {
  protected:
    T* m_data;
    size_t m_rows, m_cols;
    size_t m_stride;

    size_t index(size_t i, size_t j) const
    {
      return (ORD==RowMajor) ? i*m_stride+j : j*m_stride+i;
    }
  public:
    MatrixView() = default;
    MatrixView(const MatrixView &) = default;
        
    MatrixView (size_t rows, size_t cols, T* data)
      : m_data(data), m_rows(rows), m_cols(cols)
    {
      m_stride = (ORD==RowMajor) ? m_cols : m_rows;
    }
    
    MatrixView (size_t rows, size_t cols, size_t stride, T* data)
      : m_data(data), m_rows(rows), m_cols(cols), m_stride(stride) { }
    
    template <typename TB>
    MatrixView (const MatrixView<TB,ORD>& m2)
      : m_data(m2.data()), m_rows(m2.rows()), m_cols(m2.cols()), m_stride(m2.stride()) { }


    MatrixView& operator= (const MatrixView& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = m2(i,j);
      return *this;
    }
    
    template <typename TB>
    MatrixView& operator= (const MatExpr<TB>& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = m2(i,j);
      return *this;
    }
        
    MatrixView& operator= (T scal)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) = scal;
      return *this;
    }
        
    T* data() const { return m_data; }
    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }
    size_t stride() const { return m_stride; }
    auto shape() const { return std::array<size_t,2>{m_rows, m_cols}; }

    T& operator()(size_t i, size_t j) 
    { 
      return m_data[index(i,j)];
    }
        
    const T& operator()(size_t i, size_t j) const 
    {
      return m_data[index(i,j)];
    }


    template <typename TB>
    MatrixView& operator+= (const MatExpr<TB>& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) += m2(i,j);
      return *this;
    }
    
    template <typename TB>
    MatrixView& operator-= (const MatExpr<TB>& m2)
    {
      for (size_t i = 0; i < m_rows; i++)
        for (size_t j = 0; j < m_cols; j++)
          (*this)(i,j) -= m2(i,j);
      return *this;
    }

    MatrixView& operator*= (T scal)
    {
      for (size_t i = 0; i < rows(); i++)
        for (size_t j = 0; j < cols(); j++)
          (*this)(i,j) *= scal;
      return *this;
    }


    auto row(size_t i) const 
    {
      if constexpr (ORD == RowMajor)
        return VectorView<T>(m_cols, m_data + i*m_stride);
      else
        return VectorView<T,size_t>(m_cols, m_stride, m_data + i);
    } 

    auto col(size_t j) const 
    {
      if constexpr (ORD == RowMajor)
        return VectorView<T,size_t>(m_rows, m_stride, m_data + j);
      else
        return VectorView<T>(m_rows, m_data + j*m_stride);
    }

    auto diag() const 
    {
      return VectorView<T,size_t>(std::min(m_rows, m_cols), m_stride+1, m_data);
    } 

    auto rowRange(size_t first, size_t next) const 
    {
      return MatrixView<T, ORD>(next - first, m_cols, m_stride, m_data+index(first,0));
    }

    auto colRange(size_t first, size_t next) const 
    {
      return MatrixView<T, ORD>(m_rows, next - first, m_stride, m_data+index(0,first));
    }
  };


  template <typename T, ORDERING ORD>
  auto trans (MatrixView<T,ORD> mat)
  {
    if constexpr (ORD==RowMajor)
      return MatrixView<T,ColMajor>(mat.cols(), mat.rows(), mat.stride(), mat.data());
    else
      return MatrixView<T,RowMajor>(mat.cols(), mat.rows(), mat.stride(), mat.data());
  }
  
  template <typename T=double, ORDERING ORD=RowMajor>
  class Matrix : public MatrixView<T,ORD>
  {
    typedef MatrixView<T,ORD> BASE;
    using BASE::m_cols;
    using BASE::m_data;
    using BASE::m_rows;

  public:
    Matrix (size_t rows, size_t cols)
      : BASE(rows, cols, new T[rows*cols]) { }
          
    Matrix (const Matrix& m2)
      : BASE(m2.rows(), m2.cols(), new T[m2.rows()*m2.cols()])
    {
      *this = m2;
    }
          
    Matrix (std::initializer_list<std::initializer_list<T>> list)
      : BASE(list.size(), list.begin()->size(), new T[list.size()*list.begin()->size()])
    {
      size_t i = 0;
      for (auto row : list)
        {
          size_t j = 0;
          for (auto val : row)
            {
              (*this)(i,j) = val;
              j++;
            }
          i++;
        }
    } 

    ~Matrix() { delete[] m_data; }

    using BASE::operator=;
    Matrix& operator= (const Matrix& m2)
    {
      if (this != &m2)
        {
          assert(m_rows==m2.m_rows && m_cols==m2.m_cols);
          BASE::operator=(m2);
        }
      return *this;
    }
                    
  };

}


#endif
