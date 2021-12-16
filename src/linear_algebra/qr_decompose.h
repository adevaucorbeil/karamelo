#pragma once

#include <givens_rotation.h>

// returns q and sets matrix to r
template<typename T, typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, 3, 3>
qr_decompose(Matrix<T, 3, 3> &matrix)
{
  Matrix<V, 3, 3> q, g;

  q = givens_rotation(matrix, 3, 1);
  q = givens_rotation(q*matrix, 2, 1)*q;
  q = givens_rotation(q*matrix, 3, 2)*q;

  return q;
}