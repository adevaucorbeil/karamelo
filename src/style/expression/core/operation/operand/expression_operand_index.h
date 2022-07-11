#pragma once

#include <expression_operand.h>

class ExpressionOperandIndex:
  public ExpressionOperand<ExpressionOperandIndex>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return i;
  }
};
