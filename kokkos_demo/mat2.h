#pragma once

#include <vec2.h>

struct        mat2 {
  double      m00;
  double      m01;
  double      m10;
  double      m11;
  
  KOKKOS_INLINE_FUNCTION
              mat2():
                mat2(1, 0, 0, 1)                    {}
  KOKKOS_INLINE_FUNCTION
              mat2(double m00,
                   double m01,
                   double m10,
                   double m11):
                m00(m00),
                m01(m01),
                m10(m10),
                m11(m11)                            {}
  KOKKOS_INLINE_FUNCTION
              mat2(const mat2 &)                    = default;
  KOKKOS_INLINE_FUNCTION
  mat2 &      operator=(const mat2 &)               = default;
  KOKKOS_INLINE_FUNCTION
              mat2(mat2 &&)                         = default;
  KOKKOS_INLINE_FUNCTION
  mat2 &      operator=(mat2 &&)                    = default;

  KOKKOS_INLINE_FUNCTION
  mat2 &      operator*=(double s)                  { m00 *= s; m01 *= s; m10 *= s; m11 *= s; return *this; }
  KOKKOS_INLINE_FUNCTION
  mat2        operator* (double s)            const { return mat2(*this) *= s; }

  KOKKOS_INLINE_FUNCTION
  mat2 &      operator/=(double s)                  { m00 /= s; m01 /= s; m10 /= s; m11 /= s; return *this; }
  KOKKOS_INLINE_FUNCTION
  mat2        operator/ (double s)            const { return mat2(*this) /= s; }

  KOKKOS_INLINE_FUNCTION
  mat2 &      operator+=(const mat2 &m)             { m00 += m.m00; m01 += m.m01; m10 += m.m10; m11 += m.m11; return *this; }
  KOKKOS_INLINE_FUNCTION
  mat2        operator+ (const mat2 &m)       const { return mat2(*this) += m; }

  KOKKOS_INLINE_FUNCTION
  mat2 &      operator-=(const mat2 &m)             { m00 -= m.m00; m01 -= m.m01; m10 -= m.m10; m11 -= m.m11; return *this; }
  KOKKOS_INLINE_FUNCTION
  mat2        operator- (const mat2 &m)       const { return mat2(*this) -= m; }

  KOKKOS_INLINE_FUNCTION
  mat2        operator*(const mat2 &m)        const { return mat2(m00*m.m00 + m01*m.m10, m00*m.m01 + m01*m.m11, m10*m.m00 + m11*m.m10, m10*m.m01 + m11*m.m11); }

  KOKKOS_INLINE_FUNCTION
  mat2        transpose()                     const { return mat2(m00, m10, m01, m11); }

  KOKKOS_INLINE_FUNCTION
  vec2        operator*(const vec2 &v)        const { return vec2(m00*v.x + m01*v.y, m10*v.x + m11*v.y); }

  KOKKOS_INLINE_FUNCTION
  double &    operator()(int i,
                         int j)                     { return i? (j? m11: m10): (j? m01: m00); }

  KOKKOS_INLINE_FUNCTION
  double      operator()(int i,
                         int j)               const { return i? (j? m11: m10): (j? m01: m00); }
};

KOKKOS_INLINE_FUNCTION
mat2          operator*(double s,
                        const mat2 &m)              { return m*s; }