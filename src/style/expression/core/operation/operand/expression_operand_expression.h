#pragma once

#include <expression_operand.h>

class ExpressionOperandExpression:
  public ExpressionOperand<ExpressionOperandExpression>
{
  string name;
  Kokkos::View<double**> registers;

public:
  ExpressionOperandExpression(const string &name):
    name(name)
  {}

  void
  set_solid(const Solid &solid) override
  {
    const map<string, Expression>::iterator &it = solid.input->expressions.find(name);
    if (it == solid.input->expressions.end())
      solid.error->all(FLERR, name + " is not an expression.\n");
    registers = solid.input->expressions[name].registers;
  }

  void
  set_grid(const Grid &grid) override
  {
    const map<string, Expression>::iterator &it = grid.input->expressions.find(name);
    if (it == grid.input->expressions.end())
      grid.error->all(FLERR, name + " is not an expression.\n");
    registers = grid.input->expressions[name].registers;
  }

  KOKKOS_INLINE_FUNCTION double
  evaluate(int i) const
  {
    return registers(0, registers.extent(1) > 1? i: 0);
  }
};
