#pragma once

#include <expression_function.h>

class ExpressionFunctionSinh:
  public ExpressionFunction<ExpressionFunctionSinh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::sinh(get_value(0, i));
  }
};

class ExpressionFunctionCosh:
  public ExpressionFunction<ExpressionFunctionCosh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::cosh(get_value(0, i));
  }
};

class ExpressionFunctionTanh:
  public ExpressionFunction<ExpressionFunctionTanh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::tanh(get_value(0, i));
  }
};

class ExpressionFunctionAsinh:
  public ExpressionFunction<ExpressionFunctionAsinh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::asinh(get_value(0, i));
  }
};

class ExpressionFunctionAcosh:
  public ExpressionFunction<ExpressionFunctionAcosh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::acosh(get_value(0, i));
  }
};

class ExpressionFunctionAtanh:
  public ExpressionFunction<ExpressionFunctionAtanh, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::atanh(get_value(0, i));
  }
};
