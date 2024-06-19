#pragma once

#include <eigendecompose.h>

// sets matrix to be diagonal with singular values and returns u and v
template<typename T, std::size_t M, std::size_t N,
  typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION void
singular_value_decompose(Matrix<T, M, N> &matrix, Matrix<V, M, M> &u, Matrix<V, N, N> &v)
{
  Matrix<V, N, N> eigenvalues = matrix.transpose()*matrix;

  v = eigendecompose(eigenvalues);

  Matrix<V, M, N> matrix_v = matrix*v;
  Matrix<V, N, M> inverse;
  
  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++)
    {
      if (i == j)
      {
        matrix(i, i) = Kokkos::sqrt(eigenvalues(i, i));

        inverse(i, i) = Kokkos::abs(matrix(i, i)) < 1e-10? 0: 1/matrix(i, i);
      }
      else
      {
        matrix(i, j) = 0;
      }
    }

  u = matrix_v*inverse;
}
