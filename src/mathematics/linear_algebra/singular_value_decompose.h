#pragma once

#include <eigendecompose.h>

#include <utility>

// sets matrix to be diagonal with singular values and returns u and v
template<typename T, std::size_t M, std::size_t N,
  typename V = decltype(std::declval<T>()/std::declval<T>())>
std::pair<Matrix<V, M, M>, Matrix<V, N, N>>
singular_value_decompose(Matrix<T, M, N> &matrix)
{
  Matrix<V, N, N> eigenvalues = matrix.transpose()*matrix;

  const Matrix<V, N, N> &v = eigendecompose(eigenvalues);

  Matrix<V, M, N> matrix_v = matrix*v;
  Matrix<V, N, M> inverse;
  
  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++)
    {
      if (i == j)
      {
        matrix(i, i) = std::sqrt(eigenvalues(i, i));

        inverse(i, i) = std::abs(matrix(i, i)) < 1e-10? 0: 1/matrix(i, i);
      }
      else
      {
        matrix(i, j) = 0;
      }
    }

  return std::make_pair(matrix_v*inverse, v);
}