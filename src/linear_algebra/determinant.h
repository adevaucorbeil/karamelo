#pragma once

#include <matrix.h>

template<typename T, typename V = decltype(std::declval<T>()*std::declval<T>())>
KOKKOS_INLINE_FUNCTION V
determinant(const Matrix<T, 3, 3> &matrix)
{
  return matrix(0, 0)*(matrix(1, 1)*matrix(2, 2) - matrix(2, 1)*matrix(1, 2)) -
         matrix(0, 1)*(matrix(1, 0)*matrix(2, 2) - matrix(1, 2)*matrix(2, 0)) +
         matrix(0, 2)*(matrix(1, 0)*matrix(2, 1) - matrix(1, 1)*matrix(2, 0));
}