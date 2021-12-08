#pragma once

#include "constants.h"

class Interval {
public:
  double x0;
  double x1;

  Interval() :
    x0(NAN),
    x1(NAN)
  {}

  Interval(double x0, double x1) :
    x0(x0),
    x1(x1)
  {}

  Interval(const Interval &interval) = default;

  void add(double x) {
    if (x > x1 || isnan(x1))
      x1 = x;
    if (x < x0 || isnan(x0))
      x0 = x;
  }

  bool contains(double x) const {
    return x0 - eps < x && x < x1 + eps;
  }

  bool intercepts(const Interval &interval) const {
    return interval.x0 - eps < x1 && interval.x1 + eps > x0;
  }
};