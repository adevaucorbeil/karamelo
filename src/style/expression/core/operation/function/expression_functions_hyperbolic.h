#pragma once

#include <expression_function.h>

class ExpressionFunctionSinh:
  public ExpressionFunction<ExpressionFunctionSinh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::sinh(get_value(0, i));
  }
};

class ExpressionFunctionCosh:
  public ExpressionFunction<ExpressionFunctionCosh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::cosh(get_value(0, i));
  }
};

class ExpressionFunctionTanh:
  public ExpressionFunction<ExpressionFunctionTanh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::tanh(get_value(0, i));
  }
};

class ExpressionFunctionAsinh:
  public ExpressionFunction<ExpressionFunctionAsinh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::asinh(get_value(0, i));
  }
};

class ExpressionFunctionAcosh:
  public ExpressionFunction<ExpressionFunctionAcosh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::acosh(get_value(0, i));
  }
};

class ExpressionFunctionAtanh:
  public ExpressionFunction<ExpressionFunctionAtanh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::atanh(get_value(0, i));
  }
};
