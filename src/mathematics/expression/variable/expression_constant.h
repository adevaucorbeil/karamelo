#pragma once

#include <expression.h>

template<typename T>
class ExpressionConstant: public Expression<T>
{
  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in)
  {}

public:
  ExpressionConstant(const T &constant)
  {
    value_cache = constant;
  }
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative()
  {
    return ExpressionConstant(T{});
  }
};
