#pragma once

#include <givens_rotation.h>

// returns q and sets matrix to r
template<typename T, std::size_t N,
  typename V = decltype(std::declval<T>()/std::declval<T>())>
Matrix<V, N, N>
qr_decompose(Matrix<T, N, N> &matrix)
{
  Matrix<V, N, N> q = Matrix<V, N, N>::identity();
  
  for (int i = 1; i < N; i++)
    for (int j = 0; j < i; j++)
      q = givens_rotation(q*matrix, i, j)*q;

  matrix = q*matrix;

  return q.transpose();
}