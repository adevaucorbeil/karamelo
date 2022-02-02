#pragma once

#include <expression.h>

template<typename T>
class ExpressionQuotient: public Expression<T>
{
  Expression<T> numerator;
  Expression<T> denominator;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in)
  {
    return numerator.getValue(solid, ip, grid, in)/denominator.getValue(solid, ip, grid, in);
  }

public:
  ExpressionQuotient(const Expression<T> &numerator, const Expression<T> &denominator):
    Expression(numerator, denominator),
    numerator(numerator),
    denominator(denominator)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative()
  {
    return (first.timeDerivative()*second - first*second.timeDerivative())/second/second;
  }

  template<typename T>
  friend KOKKOS_INLINE_FUNCTION ExpressionQuotient<T>
  operator/(const Expression<T> &numerator, const Expression<T> &denominator)
  {
    return ExpressionQuotient(numerator, denominator);
  }
};
