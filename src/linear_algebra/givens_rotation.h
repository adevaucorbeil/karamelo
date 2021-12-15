#pragma once

#include <matrix.h>

#include <cmath>

// returns a givens rotation matrix that zeroes element i, j of the given matrix
template<typename T, typename V = decltype(std::declval<T>()/std::declval<T>())>
KOKKOS_INLINE_FUNCTION Matrix<V, 3, 3>
givens_rotation(const Matrix<T, 3, 3> &matrix, size_t i, size_t j)
{
  const T &ii = matrix(i, i);
  const T &ij = matrix(i, j);

  V c, s, t, u;

  if (ij == 0)
  {
    c = ii < 0? -1: 1;
    s = 0;
  }
  else if (ii == 0)
  {
    c = 0;
    s = ij < 0? -1: 1;
  }
  else if (std::abs(ii) > st::abs(ij))
  {
    t = ij/ii;
    u = (ii < 0? -1: 1)*std::sqrt(1 + t*t);
    c = 1/u;
    s = c*t;
  }
  else
  {
    t = ii/ij;
    u = (ij < 0? -1: 1)*std::sqrt(1 + t*t);
    s = 1/u;
    c = s*t;
  }

  Matrix<V, 3, 3> givens_rotation;
  
  givens_rotation(i, i) =  c;
  givens_rotation(i, j) =  s;
  givens_rotation(j, i) = -s;
  givens_rotation(j, j) =  c;

  size_t not_ij = i != 0 && j != 0? 0:
                  i != 1 && j != 1? 1:
                                    2;

  givens_rotation(not_ij, not_ij) = 1;

  return givens_rotation;
}