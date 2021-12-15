#pragma once

#include <qr_decompose.h>

// diagonalizes matrix and returns matrix of eigenvectors
template<typename T, typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, 3, 3>
eigen_decompose(Matrix<T, 3, 3> &matrix)
{
  Matrix<V, 3, 3> eigenvectors;

  return eigenvectors;
}