#pragma once

#include <Kokkos_Core.hpp>

struct         vec2 {
  double      x;
  double      y;

  KOKKOS_INLINE_FUNCTION
              vec2():
                vec2(0, 0)                          {}
  KOKKOS_INLINE_FUNCTION
              vec2(double x,
                   double y):
                x(x),
                y(y)                                {}
  KOKKOS_INLINE_FUNCTION
              vec2(const vec2 &)                    = default;
  KOKKOS_INLINE_FUNCTION
  vec2 &      operator=(const vec2 &)               = default;
  KOKKOS_INLINE_FUNCTION
              vec2(vec2 &&)                         = default;
  KOKKOS_INLINE_FUNCTION
  vec2 &      operator=(vec2 &&)                    = default;

  KOKKOS_INLINE_FUNCTION
  vec2 &      operator+=(double s)                  { x += s; y += s; return *this; }
  KOKKOS_INLINE_FUNCTION
  vec2        operator+ (double s)            const { return vec2(*this) += s; }

  KOKKOS_INLINE_FUNCTION
  vec2 &      operator-=(double s)                  { x -= s; y -= s; return *this; }
  KOKKOS_INLINE_FUNCTION
  vec2        operator- (double s)            const { return vec2(*this) -= s; }

  KOKKOS_INLINE_FUNCTION
  vec2 &      operator*=(double s)                  { x *= s; y *= s; return *this; }
  KOKKOS_INLINE_FUNCTION
  vec2        operator* (double s)            const { return vec2(*this) *= s; }

  KOKKOS_INLINE_FUNCTION
  vec2 &      operator/=(double s)                  { x /= s; y /= s; return *this; }
  KOKKOS_INLINE_FUNCTION
  vec2        operator/ (double s)            const { return vec2(*this) /= s; }

  KOKKOS_INLINE_FUNCTION
  vec2 &      operator+=(const vec2 &v)             { x += v.x; y += v.y; return *this; }
  KOKKOS_INLINE_FUNCTION
  vec2        operator+ (const vec2 &v)       const { return vec2(*this) += v; }

  KOKKOS_INLINE_FUNCTION
  vec2 &      operator-=(const vec2 &v)             { x -= v.x; y -= v.y; return *this; }
  KOKKOS_INLINE_FUNCTION
  vec2        operator- (const vec2 &v)       const { return vec2(*this) -= v; }

  KOKKOS_INLINE_FUNCTION
  vec2        square()                        const { return vec2(x*x, y*y); }

  KOKKOS_INLINE_FUNCTION
  double      norm()                          const { return hypot(x, y); }
};

KOKKOS_INLINE_FUNCTION
vec2          operator+(double s,
                        const vec2 &v)              { return v + s; }

KOKKOS_INLINE_FUNCTION
vec2          operator-(double s,
                        const vec2 &v)              { return (v - s)*-1; }

KOKKOS_INLINE_FUNCTION
vec2          operator*(double s,
                        const vec2 &v)              { return v*s; }