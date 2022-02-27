#pragma once

#include <qr_decompose.h>

// diagonalizes matrix and returns matrix of eigenvectors
template<typename T, std::size_t N,
  typename V = decltype(std::declval<T>()/std::declval<T>())>
Matrix<V, N, N>
eigendecompose(Matrix<T, N, N> &matrix)
{
  Matrix<V, N, N> eigenvectors = Matrix<V, N, N>::identity();
  Matrix<V, N, N> deviation    = Matrix<V, N, N>::identity();

  for (int i = 0; i < 500 && deviation.norm() > 1e-8; i++)
  {
    deviation = matrix;

    int k = i%(N - 1) + 1;
    const T &a = matrix(k - 1, k - 1);
    const T &b = matrix(k - 1, k    );
    const T &c = matrix(k    , k    );
    const V &d = (a - c)/2;
    Matrix<V, N, N> shift;
    double denominator = std::abs(d) + std::hypot(d, b);
    if (abs(denominator) > 1e-10)
     shift = (c - sign(d)*b*b/denominator)*Matrix<V, N, N>::identity();

    Matrix<V, N, N> r = matrix - shift;
    const Matrix<V, N, N> &q = qr_decompose(r);

    matrix = r*q + shift;
    eigenvectors = eigenvectors*q;

    deviation -= matrix;
  }

  return eigenvectors;
}
