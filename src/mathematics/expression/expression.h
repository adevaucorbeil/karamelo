#pragma once

#include <ostream>

#include <Kokkos_Macros.hpp>

class Solid;
class Grid;

template<typename T/*, int PRECEDENCE*/>
class Expression
{
  bool constant;
  bool only_time_dependent;

  virtual KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) = 0;

  virtual KOKKOS_INLINE_FUNCTION void
  print(std::ostream &stream) = 0;

protected:
  T value_cache;

  KOKKOS_INLINE_FUNCTION
  Expression():
    Expression(true, true)
  {}
  KOKKOS_INLINE_FUNCTION
  Expression(bool constant, bool only_time_dependent):
    constant(constant),
    only_time_dependent(only_time_dependent)
  {}
  template<typename ...T>
  KOKKOS_INLINE_FUNCTION
  Expression(const ExpressionBase &child, T &&...children):
    Expression(children...)
  {
    constant            &= child.constant;
    only_time_dependent &= child.only_time_dependent;
  }

  virtual KOKKOS_INLINE_FUNCTION
  ~Expression()
  {}

  KOKKOS_INLINE_FUNCTION const T &
  getValue(const Solid *solid, int ip, const Grid *grid, int in, bool time_changed)
  {
    if (constant || only_time_dependent && !time_changed)
      updateValueCache(solid, ip, grid, in);

    return value_cache;
  }

public:
  KOKKOS_INLINE_FUNCTION const T &
  operator(const Solid &solid, int ip, bool time_changed)
  {
    return getValue(&solid, ip, nullptr, 0, time_changed);
  }

  KOKKOS_INLINE_FUNCTION const T &
  operator(const Grid &grid, int in, bool time_changed)
  {
    return getValue(nullptr, 0, &grid, in, time_changed);
  }

  virtual KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() = 0;

  friend KOKKOS_INLINE_FUNCTION std::ostream &
  operator<<(std::ostream &steam, const Expression<T, PRECEDENCE> &expression)
  {
    expression.print(stream);
  }
};
