#pragma once

#include <eigendecompose.h>

#include <utility>

// sets matrix to be diagonal with singular values and returns u and v
template<typename T, typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION std::pair<Matrix<V, 3, 3>, Matrix<V, 3, 3>>
singular_value_decompose(Matrix<T, 3, 3> &matrix)
{
  std::pair<Matrix<V, 3, 3>, Matrix<V, 3, 3>> uv;

  return uv;
}