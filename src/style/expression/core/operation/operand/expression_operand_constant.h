#pragma once

#include <expression_operand.h>

template <const float &CONSTANT>
class ExpressionOperandConstant:
  public ExpressionOperand<ExpressionOperandConstant<CONSTANT>>
{
public:
  KOKKOS_INLINE_FUNCTION float
  evaluate(int i) const
  {
    return CONSTANT;
  }
};
