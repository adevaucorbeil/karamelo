#pragma once

#include <expression_base.h>

template<typename T>
class Expression: public ExpressionBase
{
  virtual KOKKOS_INLINE_FUNCTION void
  updateValueCache(const Solid *solid, int ip, const Grid *grid, int in) = 0;

protected:
  T value_cache;

  using ExpressionBase::ExpressionBase;

  KOKKOS_INLINE_FUNCTION const T &
  getValue(const Solid *solid, int ip, const Grid *grid, int in)
  {
    if (cacheNeedsUpdating(solid, grid))
      updateValueCache(solid, ip, grid, in);

    return value_cache;
  }

public:
  KOKKOS_INLINE_FUNCTION const T &
  operator(const Solid &solid, int ip)
  {
    return getValue(&solid, ip, nullptr, 0);
  }

  KOKKOS_INLINE_FUNCTION const T &
  operator(const Grid &grid, int in)
  {
    return getValue(nullptr, 0, &grid, in);
  }

  virtual KOKKOS_INLINE_FUNCTION Expression<T>
  timeDerivative() = 0;
};
