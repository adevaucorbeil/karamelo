#pragma once

#include <type_traits>
#include <cassert>
#include <ostream>
#include <cstddef>

#include <Kokkos_Macros.hpp>

template<typename S>
using enable_if_scalar_t = std::enable_if_t<std::is_arithmetic<
  std::remove_reference_t<S>>::value>;

template<typename T, std::size_t M, std::size_t N>
  class Matrix
{
public:
  T elements[M][N];

  // constructors
  KOKKOS_INLINE_FUNCTION
  Matrix():
    elements{}
  {}

  template<typename ...U, typename = std::enable_if_t<sizeof...(U) == M*N>>
  KOKKOS_INLINE_FUNCTION
  Matrix(U &&...elements):
    elements{ static_cast<T>(elements)... }
  {}

  template<typename U>
  KOKKOS_INLINE_FUNCTION
  Matrix(const Matrix<U, M, N> &matrix)
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        elements[i][j] = matrix(i, j);
  }

  template<typename U>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
  operator=(const Matrix<U, M, N> &matrix)
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        elements[i][j] = matrix(i, j);

    return *this;
  }

  // accessors
  KOKKOS_INLINE_FUNCTION const T &
  operator()(std::size_t i, std::size_t j) const
  {
    assert(i < M && j < N);

    return elements[i][j];
  }

  KOKKOS_INLINE_FUNCTION T &
  operator()(std::size_t i, std::size_t j)
  {
    assert(i < M && j < N);

    return elements[i][j];
  }

  KOKKOS_INLINE_FUNCTION const T &
  operator()(std::size_t i) const
  {
    static_assert(N == 1, "Vector access requires N == 1");

    return operator()(i, 0);
  }

  KOKKOS_INLINE_FUNCTION T &
  operator()(std::size_t i)
  {
    static_assert(N == 1, "Vector access requires N == 1");

    return operator()(i, 0);
  }

  KOKKOS_INLINE_FUNCTION const T &
  operator[](std::size_t i) const
  {
    return operator()(i);
  }

  KOKKOS_INLINE_FUNCTION T &
  operator[](std::size_t i)
  {
    return operator()(i);
  }

  KOKKOS_INLINE_FUNCTION const T &
  x() const
  {
    return operator()(0);
  }

  KOKKOS_INLINE_FUNCTION T &
  x()
  {
    return operator()(0);
  }

  KOKKOS_INLINE_FUNCTION const T &
  y() const
  {
    static_assert(M > 1, "Y requires dimension >= 2");

    return operator()(1);
  }

  KOKKOS_INLINE_FUNCTION T &
  y()
  {
    static_assert(M > 1, "Y requires dimension >= 2");

    return operator()(1);
  }

  KOKKOS_INLINE_FUNCTION const T &
  z() const
  {
    static_assert(M > 1, "Z requires dimension >= 3");

    return operator()(2);
  }

  KOKKOS_INLINE_FUNCTION T &
  z()
  {
    static_assert(M > 1, "Z requires dimension >= 3");

    return operator()(2);
  }

  KOKKOS_INLINE_FUNCTION const T &
  w() const
  {
    static_assert(M > 1, "W requires dimension >= 4");

    return operator()(3);
  }

  KOKKOS_INLINE_FUNCTION T &
  w()
  {
    static_assert(M > 1, "W requires dimension >= 4");

    return operator()(3);
  }

  // addition/subtraction
  template<typename U, typename V = decltype(std::declval<T>() + std::declval<U>())>
  KOKKOS_INLINE_FUNCTION Matrix<V, M, N>
  operator+(const Matrix<U, M, N> &matrix) const
  {
    Matrix<V, M, N> sum;
    
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        sum(i, j) = elements[i][j] + matrix(i, j);

    return sum;
  }
  
  template<typename U, typename V = decltype(std::declval<T>() - std::declval<U>())>
  KOKKOS_INLINE_FUNCTION Matrix<V, M, N>
  operator-(const Matrix<U, M, N> &matrix) const
  {
    Matrix<V, M, N> difference;
    
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        difference(i, j) = elements[i][j] - matrix(i, j);

    return difference;
  }

  template<typename U>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
  operator+=(const Matrix<U, M, N> &matrix)
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        elements[i][j] += matrix(i, j);

    return *this;
  }
  
  template<typename U>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
  operator-=(const Matrix<U, M, N> &matrix)
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        elements[i][j] -= matrix(i, j);

    return *this;
  }

  // scalar multiplication/division
  template<typename S, typename = enable_if_scalar_t<S>,
    typename U = decltype(std::declval<T>()*std::declval<S>())>
  KOKKOS_INLINE_FUNCTION Matrix<U, M, N>
  operator*(S &&scalar) const
  {
    Matrix<U, M, N> product;

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        product(i, j) = elements[i][j]*scalar;

    return product;
  }

  template<typename S, typename = enable_if_scalar_t<S>,
    typename U = decltype(std::declval<T>()/std::declval<S>())>
  KOKKOS_INLINE_FUNCTION Matrix<U, M, N>
  operator/(S &&scalar) const
  {
    Matrix<U, M, N> quotient;

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        quotient(i, j) = elements[i][j]/scalar;

    return quotient;
  }

  template<typename S, typename = enable_if_scalar_t<S>>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
  operator*=(const S &scalar)
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        elements[i][j] *= scalar;

    return *this;
  }

  template<typename S, typename = enable_if_scalar_t<S>>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
  operator/=(const S &scalar)
  {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        elements[i][j] /= scalar;

    return *this;
  }

  KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
  operator-() const
  {
    return (*this)*-1;
  }

  // matrix product
  template<typename U, std::size_t O>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, O>
  operator*(const Matrix<U, N, O> &matrix) const
  {
    Matrix<T, M, O> product;

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        for (int k = 0; k < O; k++)
          product(i, k) += elements[i][j]*matrix(j, k);

    return product;
  }

  // dot product
  template<typename U, typename V = decltype(std::declval<T>()*std::declval<U>())>
  KOKKOS_INLINE_FUNCTION V
  dot(const Matrix<U, M, N> &matrix) const
  {
    V dot_product{};

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        dot_product += elements[i][j]*matrix(i, j);

    return dot_product;
  }

  // cross product
  template<typename U, typename V = decltype(std::declval<T>()*std::declval<U>())>
  KOKKOS_INLINE_FUNCTION Matrix<V, 3, 1>
  cross(const Matrix<U, 3, 1> &vector) const
  {
    static_assert(M == 3 && N == 1, "Cross product requires vector of dimension 3");

    return Matrix<V, 3, 1>(y()*vector.z() - z()*vector.y(),
                           z()*vector.x() - x()*vector.z(),
                           x()*vector.y() - y()*vector.x());
  }

  // norms
  KOKKOS_INLINE_FUNCTION decltype(std::declval<T>()*std::declval<T>())
  norm2() const
  {
    return dot(*this);
  }

  KOKKOS_INLINE_FUNCTION decltype(std::sqrt(std::declval<T>()*std::declval<T>()))
  norm() const
  {
    return std::sqrt(norm2());
  }

  // normalize columns
  template<typename U = decltype(std::declval<T>()*std::declval<T>())>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
  unit() const
  {
    Matrix<U, M, N> unit;

    for (int i = 0; i < N; i++)
    {
      U norm{};

      for (int j = 0; j < M; j++)
        norm += elements[j][i]*elements[j][i];

      norm = std::sqrt(norm);

      for (int j = 0; j < M; j++)
        unit(j, i) = elements[j][i]/norm;
    }

    return unit;
  }
  
  template<typename U = decltype(std::declval<T>()*std::declval<T>())>
  KOKKOS_INLINE_FUNCTION Matrix<T, M, N> &
  normalize()
  {
    for (int i = 0; i < N; i++)
    {
      U norm{};

      for (int j = 0; j < M; j++)
        norm += elements[j][i]*elements[j][i];

      norm = std::sqrt(norm);

      for (int j = 0; j < M; j++)
        elements[j][i] /= norm;
    }

    return *this;
  }

  // element-wise square
  template<typename U = decltype(std::declval<T>()*std::declval<T>())>
  KOKKOS_INLINE_FUNCTION Matrix<U, M, N>
  square() const
  {
    Matrix<U, M, N> square;
    
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        square(i, j) = elements[i][j]*elements[i][j];

    return square;
  }

  // transpose
  KOKKOS_INLINE_FUNCTION Matrix<T, N, M>
  transpose() const
  {
    Matrix<T, N, M> result;

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        result(j, i) = elements[i][j];

    return result;
  }

  // identity and trace
  static KOKKOS_INLINE_FUNCTION Matrix<T, M, N>
  identity()
  {
    static_assert(M == N, "Identity requires matrix to be square");

    Matrix<T, M, N> identity;

    for (int i = 0; i < M && i < N; i++)
      identity(i, i) = 1;

    return identity;
  }

  template<typename U = decltype(std::declval<T>() + std::declval<T>())>
  KOKKOS_INLINE_FUNCTION U
  trace() const
  {
    static_assert(M == N, "Trace requires matrix to be square");

    U trace{};

    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        trace += elements[i][j];

    return trace;
  }
};

// scale
template<typename S, typename = enable_if_scalar_t<S>,
  typename T, std::size_t M, std::size_t N,
  typename U = decltype(std::declval<T>()/std::declval<S>())>
KOKKOS_INLINE_FUNCTION Matrix<U, M, N>
operator*(const S &scalar, const Matrix<T, M, N> &matrix)
{
  return matrix*scalar;
}

// printing
template<typename T, std::size_t M, std::size_t N>
KOKKOS_INLINE_FUNCTION std::ostream &
operator<<(std::ostream &os, const Matrix<T, M, N> &matrix)
{
    os << "{ ";

    for (int i = 0; i < M; i++)
    {
      os << "{ ";
      
      for (int j = 0; j < N; j++)
        os << matrix(i, j) << (j + 1 == N? "": ", ");

      os << " }" << (i + 1 == M? "": ", ");
    }

    return os << " }";
}

// typedefs
template<typename T, std::size_t N>
using Vector = Matrix<T, N, 1>;

using Vector2d = Vector<double, 2>;
using Vector3d = Vector<double, 3>;
using Vector4d = Vector<double, 4>;
using Matrix2d = Matrix<double, 2, 2>;
using Matrix3d = Matrix<double, 3, 3>;
using Matrix4d = Matrix<double, 4, 4>;