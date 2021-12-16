#pragma once

#include <qr_decompose.h>

// diagonalizes matrix and returns matrix of eigenvectors
template<typename T, typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, 3, 3>
eigendecompose(Matrix<T, 3, 3> &matrix)
{
  Matrix<V, 3, 3> eigenvectors = Matrix<V, 3, 3>::identity();

  for (int i = 0; i < 20; i++)
  {
    const Matrix<V, 3, 3> &q = qr_decompose(matrix);
    eigenvectors = q*eigenvectors;
    matrix = matrix*q;
  }

  return eigenvectors;
}