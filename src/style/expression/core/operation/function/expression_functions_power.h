#pragma once

#include <expression_function.h>

class ExpressionFunctionPow:
  public ExpressionFunction<ExpressionFunctionPow, 2>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::pow(get_value(0, i), get_value(1, i));
  }
};

class ExpressionFunctionSqrt:
  public ExpressionFunction<ExpressionFunctionSqrt, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::sqrt(get_value(0, i));
  }
};

class ExpressionFunctionCbrt:
  public ExpressionFunction<ExpressionFunctionCbrt, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::cbrt(get_value(0, i));
  }
};

class ExpressionFunctionHypot:
  public ExpressionFunction<ExpressionFunctionHypot, 2>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::hypot(get_value(0, i), get_value(1, i));
  }
};

