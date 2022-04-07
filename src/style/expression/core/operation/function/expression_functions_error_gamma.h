#pragma once

#include <expression_function.h>

class ExpressionFunctionErf:
  public ExpressionFunction<ExpressionFunctionErf, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::erf(get_value(0, i));
  }
};

class ExpressionFunctionErfc:
  public ExpressionFunction<ExpressionFunctionErfc, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::erfc(get_value(0, i));
  }
};

class ExpressionFunctionTgamma:
  public ExpressionFunction<ExpressionFunctionTgamma, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::tgamma(get_value(0, i));
  }
};

class ExpressionFunctionLgamma:
  public ExpressionFunction<ExpressionFunctionLgamma, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return Kokkos::Experimental::lgamma(get_value(0, i));
  }
};
