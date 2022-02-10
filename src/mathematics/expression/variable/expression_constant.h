#pragma once

#include <expression.h>

#include <string>

template<typename T>
class ExpressionConstant: public Expression<T>
{
  std::string name;

  KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in)
  {}

  KOKKOS_INLINE_FUNCTION void
  print(std::ostream &stream)
  {
    stream << value_cache;
  }

public:
  ExpressionConstant(const T &constant, const std::string &name):
    name(name)
  {
    value_cache = constant;
  }
  
  KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative()
  {
    return ExpressionConstant(T{});
  }
};
