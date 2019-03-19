#ifndef LMP_MATH_SPECIAL_H
#define LMP_MATH_SPECIAL_H

#include <math.h>


namespace MathSpecial {
  // x**2, use instead of pow(x,2.0)

  static inline double square(const double &x) { return x*x; }
}

#endif
