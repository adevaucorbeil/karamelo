#pragma once
//
//#include <eigendecompose.h>
//
//template<typename T, std::size_t N,
//  typename V = decltype(std::declval<T>()*std::declval<T>())>
//KOKKOS_INLINE_FUNCTION Matrix<V, N, N>
//inverse(const Matrix<T, N, N> &matrix)
//{
//  Matrix<V, N, N> eigenvalues = matrix;
//  const Matrix<V, N, N> &eigenvectors = eigendecompose(eigenvalues);
//
//  for (int i = 0; i < N; i++)
//    eigenvalues(i, i) = 1/eigenvalues(i, i);
//
//  return eigenvectors*eigenvalues*eigenvectors.transpose();
//}

#include <determinant.h>

template<typename T, typename V = decltype(std::declval<T>()*std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, 3, 3>
inverse(const Matrix<T, 3, 3> &matrix)
{
  Matrix<V, 3, 3> inverse;

  inverse(0, 0) = (matrix(1, 1)*matrix(2, 2) - matrix(2, 1)*matrix(1, 2));
  inverse(0, 1) = (matrix(0, 2)*matrix(2, 1) - matrix(0, 1)*matrix(2, 2));
  inverse(0, 2) = (matrix(0, 1)*matrix(1, 2) - matrix(0, 2)*matrix(1, 1));
  inverse(1, 0) = (matrix(1, 2)*matrix(2, 0) - matrix(1, 0)*matrix(2, 2));
  inverse(1, 1) = (matrix(0, 0)*matrix(2, 2) - matrix(0, 2)*matrix(2, 0));
  inverse(1, 2) = (matrix(1, 0)*matrix(0, 2) - matrix(0, 0)*matrix(1, 2));
  inverse(2, 0) = (matrix(1, 0)*matrix(2, 1) - matrix(2, 0)*matrix(1, 1));
  inverse(2, 1) = (matrix(2, 0)*matrix(0, 1) - matrix(0, 0)*matrix(2, 1));
  inverse(2, 2) = (matrix(0, 0)*matrix(1, 1) - matrix(1, 0)*matrix(0, 1));

  return inverse /= determinant(matrix);
}