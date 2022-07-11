#pragma once

#include <expression_function.h>

class ExpressionFunctionSin:
  public ExpressionFunction<ExpressionFunctionSin, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::sin(get_value(0, i));
  }
};

class ExpressionFunctionCos:
  public ExpressionFunction<ExpressionFunctionCos, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::cos(get_value(0, i));
  }
};

class ExpressionFunctionTan:
  public ExpressionFunction<ExpressionFunctionTan, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::tan(get_value(0, i));
  }
};

class ExpressionFunctionAsin:
  public ExpressionFunction<ExpressionFunctionAsin, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::asin(get_value(0, i));
  }
};

class ExpressionFunctionAcos:
  public ExpressionFunction<ExpressionFunctionAcos, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::acos(get_value(0, i));
  }
};

class ExpressionFunctionAtan:
  public ExpressionFunction<ExpressionFunctionAtan, 1>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::atan(get_value(0, i));
  }
};

class ExpressionFunctionAtan2:
  public ExpressionFunction<ExpressionFunctionAtan2, 2>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return Kokkos::Experimental::atan2(get_value(0, i), get_value(1, i));
  }
};
