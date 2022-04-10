#pragma once

#include <expression_operand.h>

class ExpressionOperandExpression:
  public ExpressionOperand<ExpressionOperandExpression>
{
  const Expression &expression;

public:
  ExpressionOperandExpression(const Expression &expression):
    expression(expression)
  {}

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return expression[i];
  }
};
