#pragma once

#include <expression_operand.h>

class ExpressionOperandExpression:
  public ExpressionOperand<ExpressionOperandExpression>
{
  Expression &expression;

public:
  ExpressionOperandExpression(Expression &expression):
    expression(expression)
  {}

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return expression[i];
  }
};
