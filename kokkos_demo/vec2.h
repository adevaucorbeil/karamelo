#pragma once

#include <cmath>

struct         vec2 {
  double      x;
  double      y;

              vec2():
                vec2(0, 0)                          {}
              vec2(double x,
                   double y):
                x(x),
                y(y)                                {}
              vec2(const vec2 &)                    = default;
  vec2 &      operator=(const vec2 &)               = default;
              vec2(vec2 &&)                         = default;
  vec2 &      operator=(vec2 &&)                    = default;

  vec2 &      operator+=(double s)                  { x += s; y += s; return *this; }
  vec2        operator+ (double s)            const { return vec2(*this) += s; }

  vec2 &      operator-=(double s)                  { x -= s; y -= s; return *this; }
  vec2        operator- (double s)            const { return vec2(*this) -= s; }

  vec2 &      operator*=(double s)                  { x *= s; y *= s; return *this; }
  vec2        operator* (double s)            const { return vec2(*this) *= s; }

  vec2 &      operator/=(double s)                  { x /= s; y /= s; return *this; }
  vec2        operator/ (double s)            const { return vec2(*this) /= s; }

  vec2 &      operator+=(const vec2 &v)             { x += v.x; y += v.y; return *this; }
  vec2        operator+ (const vec2 &v)       const { return vec2(*this) += v; }

  vec2 &      operator-=(const vec2 &v)             { x -= v.x; y -= v.y; return *this; }
  vec2        operator- (const vec2 &v)       const { return vec2(*this) -= v; }

  vec2        square()                        const { return vec2(x*x, y*y); }

  double      norm()                          const { return hypot(x, y); }
};

inline vec2 operator+(double s,
                      const vec2 &v)              { return v + s; }

inline vec2 operator-(double s,
                      const vec2 &v)              { return (v - s)*-1; }

inline vec2 operator*(double s,
                      const vec2 &v)              { return v*s; }