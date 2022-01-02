#pragma once

#include <eigendecompose.h>

template<typename T, std::size_t N,
  typename V = decltype(std::declval<T>()*std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, N, N>
inverse(const Matrix<T, N, N> &matrix)
{
  Matrix<V, N, N> eigenvalues = matrix;
  const Matrix<V, N, N> &eigenvectors = eigendecompose(eigenvalues);

  for (int i = 0; i < N; i++)
    eigenvalues(i, i) = 1/eigenvalues(i, i);

  return eigenvectors*eigenvalues*eigenvectors.transpose();
}