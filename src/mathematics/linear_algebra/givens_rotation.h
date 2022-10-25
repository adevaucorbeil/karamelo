#pragma once

#include <matrix.h>
#include <math_utility.h>

// returns a givens rotation matrix that zeroes element i, j of the given matrix
template<typename T, std::size_t M, std::size_t N,
  typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, M, M>
givens_rotation(const Matrix<T, M, N> &matrix, std::size_t i, std::size_t j)
{
  const T &jj =  matrix(j, j);
  const T &ij = -matrix(i, j);

  V c, s;

  if (ij == 0)
  {
    c = sign(jj);
    s = 0;
  }
  else if (jj == 0)
  {
    c = 0;
    s = sign(ij);
  }
  else if (Kokkos::Experimental::abs(jj) > Kokkos::Experimental::abs(ij))
  {
    V t = ij/jj;
    c = sign(jj)/Kokkos::Experimental::sqrt(1 + t*t);
    s = c*t;
  }
  else
  {
    V t = jj/ij;
    s = sign(ij)/Kokkos::Experimental::sqrt(1 + t*t);
    c = s*t;
  }

  Matrix<V, M, M> givens_rotation;
  
  givens_rotation(i, i) =  c;
  givens_rotation(i, j) =  s;
  givens_rotation(j, i) = -s;
  givens_rotation(j, j) =  c;

  for (int k = 0; k < M; k++)
    if (k != i && k != j)
      givens_rotation(k, k) = 1;

  return givens_rotation;
}