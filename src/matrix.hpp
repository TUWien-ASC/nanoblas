#ifndef FILE_MATRIX
#define FILE_MATRIX

#include "matexpr.hpp"

namespace nanoblas
{
  enum ORDERING { RowMajor, ColMajor };

  template <typename T, ORDERING ORD>
  class MatrixView : public MatExpr<MatrixView<T,ORD>>
  {
  protected:
    T * data_;
    size_t rows_, cols_;
    size_t dist_;
  public:
    MatrixView() = default;
    MatrixView(const MatrixView &) = default;
        
    MatrixView (size_t rows, size_t cols, T * data)
      : data_(data), rows_(rows), cols_(cols)
    {
      dist_ = (ORD==RowMajor) ? cols_ : rows_;
    }
    
    MatrixView (size_t rows, size_t cols, size_t dist, T * data)
      : data_(data), rows_(rows), cols_(cols), dist_(dist) { }
    
    template <typename TB, ORDERING ORD2>
    MatrixView (const MatrixView<TB,ORD2> & m2)
      : data_(m2.Data()), rows_(m2.Rows()), cols_(m2.Cols()) { }
        
    template <typename TB>
    MatrixView & operator= (const MatExpr<TB> & m2)
    {
      for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
          (*this)(i,j) = m2(i,j);
      return *this;
    }
        
    MatrixView & operator= (T scal)
    {
      for (size_t i = 0; i < rows_; i++)
        for (size_t j = 0; j < cols_; j++)
          (*this)(i,j) = scal;
      return *this;
    }
        
    T * data() const { return data_; }
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    size_t dist() const { return dist_; }
    
    T & operator()(size_t i, size_t j) 
    { 
      if constexpr (ORD == RowMajor)
        return data_[i*dist_ + j];
      else
        return data_[j*dist_ + i];
    }
        
    const T & operator()(size_t i, size_t j) const 
    { 
      if constexpr (ORD == RowMajor)
        return data_[i*dist_ + j];
      else
        return data_[j*dist_ + i];
    }
  };


  template <typename T, ORDERING ORD>
  auto Trans (MatrixView<T,ORD> mat)
  {
    if constexpr (ORD==RowMajor)
      return MatrixView<T,ColMajor>(mat.cols(), mat.rows(), mat.dist(), mat.data());
    else
      return MatrixView<T,RowMajor>(mat.cols(), mat.rows(), mat.dist(), mat.data());
  }
  
  template <typename T, ORDERING ORD=RowMajor>
  class Matrix : public MatrixView<T,ORD>
  {
    typedef MatrixView<T,ORD> BASE;
    using BASE::rows_;
    using BASE::cols_;
    T * own_data_;
  public:
    Matrix (size_t rows, size_t cols)
      : BASE(rows, cols, new T[rows*cols]), own_data_(BASE::data()) { }
          
    Matrix (const Matrix & m2)
      : BASE(m2.Rows(), m2.Cols(), new T[m2.Rows()*m2.Cols()]),
        own_data_(BASE::data())
    {
      *this = m2;
    }
          
    ~Matrix() { delete[] own_data_; }

    using BASE::operator=;
    Matrix & operator= (const Matrix & m2)
    {
      if (this != &m2)
        {
          if (rows_ != m2.Rows() || cols_ != m2.Cols())
            throw std::runtime_error("Matrix assignment: incompatible sizes");
          BASE::operator=(m2);
        }
      return *this;
    }
                    
  };

}


#endif
