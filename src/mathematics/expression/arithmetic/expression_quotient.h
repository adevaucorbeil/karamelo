#pragma once

#include <expression.h>

template<typename T>
class ExpressionProduct: public Expression<T>
{
  Expression<T> first;
  Expression<T> second;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in)
  {
    return first.getValue(solid, ip, grid, in)*second.getValue(solid, ip, grid, in);
  }

  KOKKOS_INLINE_FUNCTION void
  print(std::ostream &stream)
  {
    stream << first << "/" << second;
  }

public:
  ExpressionProduct(const Expression<T> &first, const Expression<T> &second):
    Expression(first, second),
    first(first),
    second(second)
  {}
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative()
  {
    return first*second.timeDerivative() + first.timeDerivative()*second;
  }

  template<typename T>
  friend KOKKOS_INLINE_FUNCTION ExpressionProduct<T>
  operator*(const Expression<T> &first, const Expression<T> &second)
  {
    return ExpressionProduct(first, second);
  }
};
