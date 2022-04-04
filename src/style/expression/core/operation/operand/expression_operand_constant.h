#pragma once

#include <expression_operand.h>

template <const double &CONSTANT>
class ExpressionOperandConstant:
  public ExpressionOperand<ExpressionOperandConstant<CONSTANT>>
{
public:
  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return CONSTANT;
  }
};
