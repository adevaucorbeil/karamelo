#pragma once

#include <Vec2.h>

struct                              Mat2 {
  double                            m00;
  double                            m01;
  double                            m10;
  double                            m11;

  KOKKOS_INLINE_FUNCTION            Mat2(double m00                     = 1,
                                         double m01                     = 0,
                                         double m10                     = 0,
                                         double m11                     = 1):
                                      m00(m00),
                                      m01(m01),
                                      m10(m10),
                                      m11(m11)                          {}

  KOKKOS_INLINE_FUNCTION            Mat2(const Mat2 &)                  = default;
  KOKKOS_INLINE_FUNCTION            Mat2(Mat2 &&)                       = default;

  KOKKOS_INLINE_FUNCTION Mat2 &     operator=(const Mat2 &)             = default;
  KOKKOS_INLINE_FUNCTION Mat2 &     operator=(Mat2 &&)                  = default;

  KOKKOS_INLINE_FUNCTION Mat2 &     operator*=(double s)                { m00 *= s; m01 *= s; m10 *= s; m11 *= s; return *this; }
  KOKKOS_INLINE_FUNCTION Mat2 &     operator/=(double s)                { m00 /= s; m01 /= s; m10 /= s; m11 /= s; return *this; }

  KOKKOS_INLINE_FUNCTION Mat2       operator* (double s)          const { return Mat2(*this) *= s; }
  KOKKOS_INLINE_FUNCTION Mat2       operator/ (double s)          const { return Mat2(*this) /= s; }

  KOKKOS_INLINE_FUNCTION Mat2       operator-()                   const { return *this*-1; }

  KOKKOS_INLINE_FUNCTION Mat2 &     operator+=(const Mat2 &m)           { m00 += m.m00; m01 += m.m01; m10 += m.m10; m11 += m.m11; return *this; }
  KOKKOS_INLINE_FUNCTION Mat2 &     operator-=(const Mat2 &m)           { m00 -= m.m00; m01 -= m.m01; m10 -= m.m10; m11 -= m.m11; return *this; }

  KOKKOS_INLINE_FUNCTION Mat2       operator+ (const Mat2 &m)     const { return Mat2(*this) += m; }
  KOKKOS_INLINE_FUNCTION Mat2       operator- (const Mat2 &m)     const { return Mat2(*this) -= m; }

  KOKKOS_INLINE_FUNCTION Mat2       operator*(const Mat2 &m)      const { return Mat2(m00*m.m00 + m01*m.m10, m00*m.m01 + m01*m.m11, m10*m.m00 + m11*m.m10, m10*m.m01 + m11*m.m11); }
  KOKKOS_INLINE_FUNCTION Vec2       operator*(const Vec2 &v)      const { return Vec2(m00*v.x + m01*v.y, m10*v.x + m11*v.y); }

  KOKKOS_INLINE_FUNCTION Mat2       transpose()                   const { return Mat2(m00, m10, m01, m11); }


  KOKKOS_INLINE_FUNCTION double &   operator()(int i,
                                               int j)                   { return i? (j? m11: m10): (j? m01: m00); }
  KOKKOS_INLINE_FUNCTION double     operator()(int i,
                                               int j)             const { return i? (j? m11: m10): (j? m01: m00); }
};

KOKKOS_INLINE_FUNCTION Mat2          operator*(double s,
                                               const Mat2 &m)           { return m*s; }