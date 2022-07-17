#pragma once

#include <expression_operator.h>

class ExpressionSum:
  public ExpressionOperatorBinary<ExpressionSum, 6>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return get_value(0, i) + get_value(1, i);
  }
};

class ExpressionDifference:
  public ExpressionOperatorBinary<ExpressionDifference, 6>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return get_value(0, i) - get_value(1, i);
  }
};

class ExpressionProduct:
  public ExpressionOperatorBinary<ExpressionProduct, 5>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return get_value(0, i)*get_value(1, i);
  }
};

class ExpressionQuotient:
  public ExpressionOperatorBinary<ExpressionQuotient, 5>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return get_value(0, i)/get_value(1, i);
  }
};

class ExpressionNegation:
  public ExpressionOperation<ExpressionNegation, false, false, 3, 1>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return -get_value(0, i);
  }
};