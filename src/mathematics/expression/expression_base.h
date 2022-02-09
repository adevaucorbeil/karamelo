#pragma once

#include <mpmtype.h>

#include <Kokkos_Macros.hpp>

class Solid;
class Grid;

class ExpressionBase
{
  bool constant;
  bool constant_if_constant_time;

  bigint t_cache = 0;

  const Solid *solid_cache = nullptr;
  int ip_cache = 0;
  const Grid *grid_cache = nullptr;
  int in_cache = 0;

  KOKKOS_INLINE_FUNCTION bool
  cacheNeedsUpdating(const Solid *solid, int ip, const Grid *grid, int in);

protected:
  KOKKOS_INLINE_FUNCTION
  ExpressionBase():
    ExpressionBase(true, true)
  {}
  KOKKOS_INLINE_FUNCTION
  ExpressionBase(bool constant, bool constant_if_constant_time):
    constant(constant),
    constant_if_constant_time(constant_if_constant_time)
  {}
  template<typename ...T>
  KOKKOS_INLINE_FUNCTION
  ExpressionBase(const ExpressionBase &child, T &&...children):
    ExpressionBase(children...)
  {
    constant                  &= child.constant;
    constant_if_constant_time &= child.constant_if_constant_time;
  }

  virtual KOKKOS_INLINE_FUNCTION
  ~ExpressionBase()
  {}
};
