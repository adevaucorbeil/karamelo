#pragma once

#include <constants.h>

#include <cmath>

#include <cmath>

using namespace std;

class Interval {
public:
  float x0;
  float x1;

  Interval() :
    x0(NAN),
    x1(NAN)
  {}

  Interval(float x0, float x1) :
    x0(x0),
    x1(x1)
  {}

  Interval(const Interval &interval) = default;

  void add(float x) {
    if (x > x1 || std::isnan(x1))
      x1 = x;
    if (x < x0 || std::isnan(x0))
      x0 = x;
  }

  bool contains(float x) const {
    return x0 - eps < x && x < x1 + eps;
  }

  bool intercepts(const Interval &interval) const {
    return interval.x0 - eps < x1 && interval.x1 + eps > x0;
  }
};
