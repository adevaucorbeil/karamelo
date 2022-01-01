#pragma once

#include <qr_decompose.h>

// diagonalizes matrix and returns matrix of eigenvectors
template<typename T, size_t N,
  typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, N, N>
eigendecompose(Matrix<T, N, N> &matrix)
{
  Matrix<V, N, N> eigenvectors = Matrix<V, N, N>::identity();
  Matrix<V, N, N> deviation    = Matrix<V, N, N>::identity();

  while (deviation.norm() > 1e-10)
  {
    deviation = matrix;

    const Matrix<V, N, N> &q = qr_decompose(matrix);
    matrix = matrix*q;

    eigenvectors = eigenvectors*q;

    deviation -= matrix;
  }

  return eigenvectors;
}