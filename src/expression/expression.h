#pragma once

#include <expression_calculator.h>

class Solid;
class Grid;

class Expression
{
  bool constant;
  bool only_time_dependent;

  double value_cache;

  /*
  KOKKOS_INLINE_FUNCTION
  Expression():
    Expression(true, true)
  {}
  KOKKOS_INLINE_FUNCTION
  Expression(bool constant, bool only_time_dependent):
    constant(constant),
    only_time_dependent(only_time_dependent)
  {}
  template<typename ...U>
  KOKKOS_INLINE_FUNCTION
  Expression(const ExpressionBase &child, U &&...children):
    Expression(children...)
  {
    constant            &= child.constant;
    only_time_dependent &= child.only_time_dependent;
  }
  */

  KOKKOS_INLINE_FUNCTION const T &
  getValue(const Solid *solid, int ip, const Grid *grid, int in, bool time_changed) const
  {
    if (!constant && (!only_time_dependent || time_changed)
      value_cache = ExpressionCalculator::calculate()

    return value_cache;
  }

public:
  KOKKOS_INLINE_FUNCTION const T &
  operator()(const Solid &solid, int ip, bool time_changed) const
  {
    return getValue(&solid, ip, nullptr, 0, time_changed);
  }

  KOKKOS_INLINE_FUNCTION const T &
  operator()(const Grid &grid, int in, bool time_changed) const
  {
    return getValue(nullptr, 0, &grid, in, time_changed);
  }
};
