#pragma once

#include <expression.h>

/*template<typename T>
class ExpressionSin: public Expression<T>
{
  Expression<T> angle;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in)
  {
    return sin(angle.getValue(solid, ip, grid, in));
  }

  KOKKOS_INLINE_FUNCTION void
  print(std::ostream &stream)
  {
    stream << "sin(" << angle << ")";
  }

public:
  ExpressionSum(const Expression<T> &angle):
    Expression(angle),
    angle(angle)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative()
  {
    return ExpressionCos(angle);
  }
};*/
