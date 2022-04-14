#pragma once

#include <expression_operand.h>

class ExpressionOperandExpression:
  public ExpressionOperand<ExpressionOperandExpression>
{
  Kokkos::View<double**> registers;

public:
  ExpressionOperandExpression(const Expression &expression):
    registers(expression.registers)
  {}

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return registers(0, registers.extent(1) > 1? i: 0);
  }
};
