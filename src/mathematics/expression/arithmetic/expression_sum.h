#pragma once

#include <expression.h>

template<typename T>
class ExpressionSum: public Expression<T>
{
  Expression<T> first;
  Expression<T> second;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in)
  {
    return first.getValue(solid, ip, grid, in) + second.getValue(solid, ip, grid, in);
  }

public:
  ExpressionSum(const Expression<T> &first, const Expression<T> &second):
    Expression(first, second),
    first(first),
    second(second)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative()
  {
    return first.timeDerivative() + second.timeDerivative();
  }

  template<typename T>
  friend KOKKOS_INLINE_FUNCTION ExpressionSum<T>
  operator+(const Expression<T> &first, const Expression<T> &second)
  {
    return ExpressionSum(first, second);
  }
};
