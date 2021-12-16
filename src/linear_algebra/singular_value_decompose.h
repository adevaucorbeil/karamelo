#pragma once

#include <eigendecompose.h>

#include <utility>

// sets matrix to be diagonal with singular values and returns u and v
template<typename T, typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION std::pair<Matrix<V, 3, 3>, Matrix<V, 3, 3>>
singular_value_decompose(Matrix<T, 3, 3> &matrix)
{
  Matrix<V, 3, 3> inner_product = matrix.transpose()*matrix;
  matrix = matrix*matrix.transpose();

  const Matrix<V, 3, 3> &u = eigendecompose(inner_product);
  const Matrix<V, 3, 3> &v = eigendecompose(matrix);
  
  for (int i = 0; i < 3; i++)
    matrix(i, i) = std::sqrt(matrix(i, i));

  return make_pair(u, v);
}