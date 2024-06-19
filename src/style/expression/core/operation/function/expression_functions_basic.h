#pragma once

#include <expression_function.h>

class ExpressionFunctionAbs:
  public ExpressionFunction<ExpressionFunctionAbs, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::abs(get_value(0, i));
  }
};

class ExpressionFunctionRemainder:
  public ExpressionFunction<ExpressionFunctionRemainder, 2>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::remainder(get_value(0, i), get_value(1, i));
  }
};

class ExpressionFunctionMax:
  public ExpressionFunction<ExpressionFunctionMax, 2>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::max(get_value(0, i), get_value(1, i));
  }
};

class ExpressionFunctionMin:
  public ExpressionFunction<ExpressionFunctionMin, 2>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::min(get_value(0, i), get_value(1, i));
  }
};
