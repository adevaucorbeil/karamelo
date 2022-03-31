#pragma once

#include <expression_operand.h>

class ExpressionOperandConstant:
  public ExpressionOperand<ExpressionOperandConstant>
{
  double value;

public:
  ExpressionOperandConstant(double value):
    value(value)
  {}

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return value;
  }
};
