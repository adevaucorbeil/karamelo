#pragma once

#include <expression_operand.h>

class ExpressionOperandLiteral:
  public ExpressionOperand<ExpressionOperandLiteral>
{
  float value;

  friend class Input;

public:
  ExpressionOperandLiteral(float value):
    value(value)
  {}

  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return value;
  }
};
