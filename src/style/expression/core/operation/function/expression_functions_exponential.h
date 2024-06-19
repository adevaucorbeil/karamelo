#pragma once

#include <expression_function.h>

class ExpressionFunctionExp:
  public ExpressionFunction<ExpressionFunctionExp, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::exp(get_value(0, i));
  }
};

class ExpressionFunctionExp2:
  public ExpressionFunction<ExpressionFunctionExp2, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::exp2(get_value(0, i));
  }
};

class ExpressionFunctionExpm1:
  public ExpressionFunction<ExpressionFunctionExpm1, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::expm1(get_value(0, i));
  }
};

class ExpressionFunctionLog:
  public ExpressionFunction<ExpressionFunctionLog, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::log(get_value(0, i));
  }
};

class ExpressionFunctionLog10:
  public ExpressionFunction<ExpressionFunctionLog10, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::log10(get_value(0, i));
  }
};

class ExpressionFunctionLog2:
  public ExpressionFunction<ExpressionFunctionLog2, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::log2(get_value(0, i));
  }
};

class ExpressionFunctionLog1p:
  public ExpressionFunction<ExpressionFunctionLog1p, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::log1p(get_value(0, i));
  }
};
