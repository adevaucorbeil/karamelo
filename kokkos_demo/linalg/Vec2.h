#pragma once

#include <Kokkos_Core.hpp>

struct                              Vec2 {
  double                            x;
  double                            y;

  KOKKOS_INLINE_FUNCTION            Vec2(double x                       = 0,    
                                         double y                       = 0):   
                                    x(x),   
                                    y(y)                                {}

  KOKKOS_INLINE_FUNCTION            Vec2(const Vec2 &)                  = default;
  KOKKOS_INLINE_FUNCTION            Vec2(      Vec2 &&)                 = default;

  KOKKOS_INLINE_FUNCTION Vec2 &     operator=(const Vec2 &)             = default;
  KOKKOS_INLINE_FUNCTION Vec2 &     operator=(      Vec2 &&)            = default;

  KOKKOS_INLINE_FUNCTION Vec2 &     operator+=(double s)                { x += s; y += s; return *this; }
  KOKKOS_INLINE_FUNCTION Vec2 &     operator-=(double s)                { x -= s; y -= s; return *this; }
  KOKKOS_INLINE_FUNCTION Vec2 &     operator*=(double s)                { x *= s; y *= s; return *this; }
  KOKKOS_INLINE_FUNCTION Vec2 &     operator/=(double s)                { x /= s; y /= s; return *this; }

  KOKKOS_INLINE_FUNCTION Vec2       operator+ (double s)          const { return Vec2(*this) += s; }
  KOKKOS_INLINE_FUNCTION Vec2       operator- (double s)          const { return Vec2(*this) -= s; }
  KOKKOS_INLINE_FUNCTION Vec2       operator/ (double s)          const { return Vec2(*this) /= s; }
  KOKKOS_INLINE_FUNCTION Vec2       operator* (double s)          const { return Vec2(*this) *= s; }

  KOKKOS_INLINE_FUNCTION Vec2 &     operator+=(const Vec2 &v)           { x += v.x; y += v.y; return *this; }
  KOKKOS_INLINE_FUNCTION Vec2       operator+ (const Vec2 &v)     const { return Vec2(*this) += v; }

  KOKKOS_INLINE_FUNCTION Vec2 &     operator-=(const Vec2 &v)           { x -= v.x; y -= v.y; return *this; }
  KOKKOS_INLINE_FUNCTION Vec2       operator- (const Vec2 &v)     const { return Vec2(*this) -= v; }

  KOKKOS_INLINE_FUNCTION Vec2       operator-()                   const { return *this*-1; }

  KOKKOS_INLINE_FUNCTION Vec2       square()                      const { return Vec2(x*x, y*y); }

  KOKKOS_INLINE_FUNCTION double     norm()                        const { return hypot(x, y); }
};

KOKKOS_INLINE_FUNCTION Vec2         operator+(double s,   
                                              const Vec2 &v)            { return v + s; }
KOKKOS_INLINE_FUNCTION Vec2         operator-(double s,   
                                              const Vec2 &v)            { return -(v - s); }
KOKKOS_INLINE_FUNCTION Vec2         operator*(double s,   
                                              const Vec2 &v)            { return v*s; }
