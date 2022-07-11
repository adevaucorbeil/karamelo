#pragma once

#include <expression_function.h>

class ExpressionFunctionCeil:
  public ExpressionFunction<ExpressionFunctionCeil, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::ceil(get_value(0, i));
  }
};

class ExpressionFunctionFloor:
  public ExpressionFunction<ExpressionFunctionFloor, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::floor(get_value(0, i));
  }
};

class ExpressionFunctionTrunc:
  public ExpressionFunction<ExpressionFunctionTrunc, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::trunc(get_value(0, i));
  }
};
