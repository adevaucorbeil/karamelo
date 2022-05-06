#pragma once

#include <expression_operand.h>

class ExpressionOperandLiteral:
  public ExpressionOperand<ExpressionOperandLiteral>
{
  double value;

public:
  ExpressionOperandLiteral(double value):
    value(value)
  {}

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return value;
  }
};