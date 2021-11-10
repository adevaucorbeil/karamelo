#pragma once

struct          Vec2 {
  double        x;
  double        y;

                Vec2(double x                       = 0,    
                     double y                       = 0):   
                  x(x),   
                  y(y)                              {}

                Vec2(const Vec2 &)                  = default;
                Vec2(      Vec2 &&)                 = default;

  Vec2 &        operator=(const Vec2 &)             = default;
  Vec2 &        operator=(      Vec2 &&)            = default;

  Vec2 &        operator+=(double s)                { x += s; y += s; return *this; }
  Vec2 &        operator-=(double s)                { x -= s; y -= s; return *this; }
  Vec2 &        operator*=(double s)                { x *= s; y *= s; return *this; }
  Vec2 &        operator/=(double s)                { x /= s; y /= s; return *this; }

  Vec2          operator+ (double s)          const { return Vec2(*this) += s; }
  Vec2          operator- (double s)          const { return Vec2(*this) -= s; }
  Vec2          operator/ (double s)          const { return Vec2(*this) /= s; }
  Vec2          operator* (double s)          const { return Vec2(*this) *= s; }

  Vec2 &        operator+=(const Vec2 &v)           { x += v.x; y += v.y; return *this; }
  Vec2          operator+ (const Vec2 &v)     const { return Vec2(*this) += v; }

  Vec2 &        operator-=(const Vec2 &v)           { x -= v.x; y -= v.y; return *this; }
  Vec2          operator- (const Vec2 &v)     const { return Vec2(*this) -= v; }

  Vec2          operator-()                   const { return *this*-1; }

  Vec2          square()                      const { return Vec2(x*x, y*y); }

  double        norm()                        const { return hypot(x, y); }
};

Vec2            operator+(double s,   
                          const Vec2 &v)            { return v + s; }
Vec2            operator-(double s,   
                          const Vec2 &v)            { return -(v - s); }
Vec2            operator*(double s,   
                          const Vec2 &v)            { return v*s; }